#include "general_libraries.h"
#include "cell_abm.h"
#include "libmesh/fe_interface.h"
//###########################################################################
//                    0 - Dead Cells
//                    1 - Tumor Cells
//                    2 - Proliferative Tumor Cells
//                    3 - Hypoxic Tumor Cells
//                    4 - Dying Tumor Cells
//                    5 - G1 Tumor Cells
//                    6 - Normoxic Cells
//###########################################################################

double max(double A, double B)
{
  if(A>B) return A;
  else return B;  
}

double point_distance(const Point& A,const Point& B)
{
  return sqrt(pow(A(0)-B(0),2)+pow(A(1)-B(1),2));
}

void update_states(list<Cell> & Cells_local, int time, Ran & ran, EquationSystems & equation_systems, const int outside_cells, int & total_tumor, int & num_dead)
{
  // ********** Define input parameters **********
  const Real N_radius     = equation_systems.parameters.get<Real>("nucleus_radius");
  const Real C_radius     = equation_systems.parameters.get<Real>("cell_radius");
  const Real A_radius     = equation_systems.parameters.get<Real>("action_radius");
  const Real tau_A        = equation_systems.parameters.get<Real>("apop_time");
  const Real tau_P        = equation_systems.parameters.get<Real>("cellc_time");
  const Real tau_g1       = equation_systems.parameters.get<Real>("g1_time");
  const Real delta_tt     = equation_systems.parameters.get<Real>("delta_tt");
  const Real alfa_p_barra = equation_systems.parameters.get<Real>("prol_intes");
  const Real alfa_a_barra = equation_systems.parameters.get<Real>("apop_intes");
  //const Real tau_NL       = equation_systems.parameters.get<Real>("lysing_time");
  //const Real f_NS         = equation_systems.parameters.get<Real>("f_NS");
  const Real sigma_h      = equation_systems.parameters.get<Real>("hypoxic_thrs");
  const double height     = equation_systems.parameters.get<Real>("domain_diameter");
  const Real gamma_A      = equation_systems.parameters.get<Real>("gamma_A");
  const Real gamma_P      = equation_systems.parameters.get<Real>("gamma_P");
  
  //********** Bookkeeping information for cells and the mesh **********
  //std::vector<Number> soln;  
  //equation_systems.build_solution_vector(soln);
  const MeshBase& mesh = equation_systems.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & system = equation_systems.get_system<TransientLinearImplicitSystem>("Diffusion");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  std::vector<dof_id_type> dof_indices;
  
  //********** Perform updates on all cells **********
  std::list<Cell>::iterator it;
  it = Cells_local.begin();
  while(it != Cells_local.end()){
    double tau = time - (*it).time;
    // ********** Tumor cells **********
    if((*it).state == 1){
      // ********** Compute nutrient concentration **********
      double sigma_cel = 0.0;
	  Point p((*it).x, (*it).y, 0.);
	  // Get point locator structure
	  UniquePtr<PointLocatorBase> point_locator = mesh.sub_point_locator();
	  // Get element where point p is
	  const Elem* top_elem = (*point_locator)(p);
	  if(top_elem){
	    // Get the location (on the reference element) of the point p located in physical space
	    Point p_master = FEInterface::inverse_map(dim, fe_type, top_elem, p);
	    dof_map.dof_indices (top_elem, dof_indices);
	    fe->reinit (top_elem);
	    const unsigned int dof_size = dof_indices.size();
	    std::vector<Real> point_phi(dof_size);
	    // Get shape function
	    for (unsigned int i=0; i != dof_size; i++)
	      point_phi[i] = FEInterface::shape(dim, fe_type, top_elem, i, p_master);
	    // Compute solution at point p
	    sigma_cel = 0.0;
	    for(unsigned int l=0; l<dof_size; l++)
	      sigma_cel += point_phi[l]*system.old_solution(dof_indices[l]);
	  }
      // Base parameter calculation
      double alfa_p = max(alfa_p_barra * ( (sigma_cel - sigma_h) / (1.0 - sigma_h)),0.0);
      double alfa_a = max(alfa_a_barra + gamma_A/(1.0+exp(-2.0*gamma_P*(sigma_h-sigma_cel))),0.0);
	  double prob_mit = (1.0-exp(-delta_tt*alfa_p))/2.0;
	  double prob_apo = (1.0-exp(-delta_tt*alfa_a))/2.0;
	  double rand_num = ran.doub();
	  if(rand_num >= 0.5){
	    // ********** Mitosis **********
	    if(prob_mit >= ran.doub()){
	      (*it).time = time;
	      (*it).prev_state = (*it).state;
	      (*it).state = 2;
	    }
	    // ********** Apoptosis **********
	    else if(prob_apo >= ran.doub()){
	      (*it).time = time;
	      (*it).prev_state = (*it).state;
	      (*it).state = 4;
	    }
	  }
	  else{
	    // ********** Apoptosis **********
	    if(prob_apo >= ran.doub()){
	      (*it).time = time;
	      (*it).prev_state = (*it).state;
	      (*it).state = 4;
	    }
	    // ********** Mitosis **********
	    else if(prob_mit >= ran.doub()){
	      (*it).time = time;
	      (*it).prev_state = (*it).state;
	      (*it).state = 2;
	    }
	  }
	  it++;
	}
    // ********** Proliferative cells **********
    else if((*it).state == 2){
	  if(tau >= (tau_P-tau_g1))
        divide_cell(Cells_local,it,ran,height);
	  it++;
	}
    // ********** Cell growth **********
    else if((*it).state == 5){
	  if(tau < tau_P){
        (*it).N_radius = (N_radius/sqrt(2.0))*sqrt(1+ (tau_g1+(tau -tau_P))/tau_g1);
	    (*it).C_radius = (C_radius/sqrt(2.0))*sqrt(1+ (tau_g1+(tau -tau_P))/tau_g1);
	    (*it).A_radius = (*it).C_radius*1.214;
	  }
	  else{
	    (*it).prev_state = 2;
	    (*it).N_radius = N_radius;
	    (*it).C_radius = C_radius;
	    (*it).A_radius = A_radius;
	    (*it).state = 1;
	    (*it).time = time;
	  }
	  it++;
	}
    // ********** Update apoptotic size **********
    else if((*it).state == 4){
	  if(tau >= tau_A)
	    it = Cells_local.erase(it);
	  else{
	    it++;
	    // ********** !new! Reduce dying cells size **********
	    //(*it).C_radius -= (*it).C_radius/(tau_A-tau+1);
	    //(*it).N_radius -= (*it).N_radius/(tau_A-tau+1);
	  }
	}
  }
}

// ********** Calculate normal **********
void normal(double& N_x,double& N_y,double x_cel,double y_cel,double height){
  // ?How does it work? !No one knows!
  if (x_cel < height/2 && x_cel != 0)
    {
      //coeficientes da reta que passa pela célula e o centro da meia circunferência        
      double a = (height/2 - y_cel)/(height/2 - x_cel);
      double b = y_cel - a*x_cel;
      //calcula o ponto que intersepta a meia circunferência com a reta acima (origem da normal)
      //Obs: sinal da raiz é negativo pois o interesse é no menor x (a esquerda do ponto central)
      double p_x = (-2*a*b+a*height+height - sqrt(-4*a*b*height+2*a*height*height-4*b*b+4*b*height))/(2*a*a+2);
      double p_y = a*p_x + b;
      //calcular vetor normal vezes distância
      N_x = x_cel - p_x;
      N_y = y_cel - p_y;
    }
  
  if (x_cel > height/2)
    {
      //coeficientes da reta que passa pela célula e o centro da meia circunferência        
      double a = (height/2 - y_cel)/(height/2 - x_cel);
      double b = y_cel - a*x_cel;
      //calcula o ponto que intersepta a meia circunferência com a reta acima (origem da normal)
      //Obs: sinal da raiz é positivo pois o interesse é no maior x (a direita do ponto central)
      double p_x = (-2*a*b+a*height+height + sqrt(-4*a*b*height+2*a*height*height-4*b*b+4*b*height))/(2*a*a+2);
      double p_y = a*p_x + b;
      //calcular vetor normal vezes distância
      N_x = x_cel - p_x;
      N_y = y_cel - p_y;
    }
  
  if (x_cel == height/2)
    {
      //calcular vetor normal vezes distância
      N_x = 0;
      if (y_cel > height/2)  N_y = y_cel-height;
      else   N_y = y_cel;
    }
}

void compute_forces(list<Cell>& Cells_local, EquationSystems& equation_systems, double height, int& outside_cells, int& total_tumor)
{
  // ********** Parameters **********
  const double c_ccr             = equation_systems.parameters.get<Real>("c_ccr");
  const double c_tta             = equation_systems.parameters.get<Real>("c_tta");
  const double c_hha             = equation_systems.parameters.get<Real>("c_hha");
  const double max_c_radius      = equation_systems.parameters.get<Real>("cell_radius");
  const double h_bin             = 3.0*max_c_radius;
  const unsigned int number_bins = ceil(height/h_bin);
  const unsigned int total_bins  = number_bins*number_bins;
  double iter                    = 0.0;
  const double delta_tt_max      = 60.0;
  // ********** Loop **********
  while (iter < delta_tt_max){
    int n =1; int m = 1; double M = 1;
    double phi_x, phi_y, psi_x, psi_y;
    // ********** Generate bins **********
    vector< list <Cell *>> Cell_Bins(total_bins);
    std::list<Cell>::iterator it;
    it = Cells_local.begin();
    while(it != Cells_local.end()){
	  unsigned int ix = floor((*it).x/h_bin);
	  unsigned int jy = floor((*it).y/h_bin);
	  unsigned int xy = ix+jy*number_bins;
	  if(ix<0 || jy<0 ||jy>=number_bins || ix>=number_bins || xy >=total_bins){
        cout << "Error" << endl;
	    cout << "Cell = ( " << (*it).x << " , " << (*it).y << " ) = ( " << ix << " , " << jy << " ) = " << xy << endl;
	    cout << "State = " << (*it).state << endl;
	    cout << "total_bins  = " << total_bins << endl;
	    cout << "number_bins = " << number_bins << endl;
	    cout << "h_bin       = " << h_bin << endl;
	    it = Cells_local.erase(it);
	  }
	  else{
	    Cell_Bins[xy].push_back(&*it);
	    it++;
	  }
	}
    // ********** Cell-cell force **********
    for(unsigned int bin_i = 0; bin_i < Cell_Bins.size(); bin_i++){
	  std::list<Cell*>::iterator cell_a;
	  for(cell_a = Cell_Bins[bin_i].begin(); cell_a != Cell_Bins[bin_i].end(); ++cell_a){
        unsigned int jy = floor(bin_i/number_bins);
	    unsigned int ix = bin_i%number_bins;
	    unsigned int bin_j;
	    vector<int> xxx{-1,0,0,1,1};
	    vector<int> yyy{1,0,1,0,1};
	    for(unsigned int p=0; p<xxx.size(); p++){
	      if(ix+xxx[p]<number_bins && ix+xxx[p]>=0 && jy+yyy[p]<number_bins && jy+yyy[p]>=0){
	        bin_j = (ix+xxx[p])+(jy+yyy[p])*number_bins;
	        std::list<Cell *>::iterator cell_b;
	        if(bin_i == bin_j){
	          cell_b = cell_a;
	          cell_b++;
	        }
	        else
	          cell_b = Cell_Bins[bin_j].begin();
            while(cell_b != Cell_Bins[bin_j].end()){
              double F_cca_i[2] = {0.0,0.0}, F_ccr_i[2] = {0.0,0.0};
	          double r_x = (*(*cell_b)).x - (*(*cell_a)).x;
	    	  double r_y = (*(*cell_b)).y - (*(*cell_a)).y;
	    	  double R_A = (*(*cell_b)).A_radius + (*(*cell_a)).A_radius;
    	      double R_N = (*(*cell_b)).N_radius + (*(*cell_a)).N_radius;
	    	  double R   = (*(*cell_b)).C_radius + (*(*cell_a)).C_radius;
	    	  potential_adh(phi_x,phi_y,r_x,r_y,R_A,n);
        	  potential_rep(psi_x,psi_y,r_x,r_y,R_N,R,M,m);
	          double c_cca;
	          if((*(*cell_a)).state != 6 && (*(*cell_b)).state != 6)
	            c_cca = c_tta;
    	      else
	            c_cca = c_hha;
	          F_cca_i[0] += -c_cca*phi_x;
	          F_cca_i[1] += -c_cca*phi_y;
	          F_ccr_i[0] += -c_ccr*psi_x;
	          F_ccr_i[1] += -c_ccr*psi_y;
	          (*(*cell_a)).F[0] += (F_cca_i[0] + F_ccr_i[0]);
	          (*(*cell_a)).F[1] += (F_cca_i[1] + F_ccr_i[1]);
	          (*(*cell_b)).F[0] -= (F_cca_i[0] + F_ccr_i[0]);
	          (*(*cell_b)).F[1] -= (F_cca_i[1] + F_ccr_i[1]);
              cell_b++;
            }
	      }
	    }
	    // ********** Cell-boundary force **********
	    double F_ct[2] = {0.0, 0.0}, F_rct[2] = {0.0, 0.0};
	    double N_x, N_y; 
	    double K = (double) 10.0;
	    double c_ct = 10.0 * K;
	    double c_rct = 4.88836 * K;
	    normal(N_x, N_y, (*(*cell_a)).x, (*(*cell_a)).y, height);
	    potential_adh(phi_x, phi_y, N_x, N_y, (*(*cell_a)).A_radius, n);
	    potential_rep(psi_x, psi_y, N_x, N_y, (*(*cell_a)).N_radius, (*(*cell_a)).C_radius, M, m);
	    F_ct[0] = -c_ct*phi_x;
	    F_ct[1] = -c_ct*phi_y;
	    F_rct[0] = -c_rct*psi_x;
	    F_rct[1] = -c_rct*psi_y;
	    (*(*cell_a)).F[0] += (F_ct[0] + F_rct[0]);
	    (*(*cell_a)).F[1] += (F_ct[1] + F_rct[1]);  
	  }
	}
    // ********** Compute cell speed **********
    double v_max  = 5.0e-3;
    double v_mean = 0.;
    double v_std  = 0.;
    double * speeds = new double [Cells_local.size()];
    unsigned int k = 0;
    for(it = Cells_local.begin(); it != Cells_local.end(); it++){
	  (*it).v[0] = -0.5*((*it).F[0]);
	  (*it).v[1] = -0.5*((*it).F[1]);
	  double max_speed = euclid_norm((*it).v[0], (*it).v[1]);
	  if(max_speed >= v_max)
	    v_max = max_speed;
	  speeds[k] = max_speed;
	  k++;
	  v_mean += max_speed;
	}
    v_mean = v_mean / Cells_local.size();
    k = 0;
    for (k = 0; k < Cells_local.size(); k++)
	  v_std += (speeds[k] - v_mean) * (speeds[k] - v_mean);
    delete[] speeds;
    v_std = sqrt(v_std / (Cells_local.size() - 1));
    // If standard deviation is low enough, take larger delta_tt value
    // This is because a low standard deviation value for velocity implies cells are in equilibrium    
    double delta_tt;
    double std_thresh = 0.001;
    if (v_std < std_thresh){
	  if(v_std == 0)
	    delta_tt = delta_tt_max;
	  else
	    delta_tt = (1.0/v_max) * (1. + std_thresh / v_std);
    }
    else
	  delta_tt = (1.0/v_max);
    if(delta_tt>delta_tt_max)
	  delta_tt = delta_tt_max;
    
    vector<int> remove_cells;
    // ********** Update cell position **********
    it = Cells_local.begin();
    while(it != Cells_local.end()){
	  (*it).x += delta_tt * (*it).v[0];
	  (*it).y += delta_tt * (*it).v[1];
	  (*it).F[0] = 0.0;
	  (*it).F[1] = 0.0;
	  double dist = distance((*it), height/2.0, height/2.0);
	  if(dist >= height/2.0){
	    it = Cells_local.erase(it);
	    outside_cells++;
	  }
	  else
	    it++;
	}
    iter += delta_tt;
  }
}

// ********** Calculate Euclidean norm between two points **********
double euclid_norm(double x,double y)
{
  return sqrt(pow(x,2.0) + pow(y,2.0));
}

// ********** Adhesion potential function **********
void potential_adh(double& phi_x,double& phi_y,double r_x,double r_y,double R_A,int n)
{
  double r = euclid_norm(r_x,r_y);
  if((r >= 0) && (r <= R_A))
    {
      double var = pow((1.0 - r/R_A),n+1)/r;
      phi_x = var*r_x;
      phi_y = var*r_y;
    }
  else
    {
      phi_x = 0;
      phi_y = 0;
    }
}

// ********** Repulsion potential function **********
void potential_rep(double& psi_x,double& psi_y,double r_x,double r_y,double R_N,double R,double M,int m)
{
  double r = euclid_norm(r_x,r_y);
  if((r >= 0) && (r < R_N))
    {
      double c = pow((1.0 - R_N/R),m+1) - M;
      double var = -(c*r/R_N + M)/r;
      psi_x = var*r_x;
      psi_y = var*r_y;
    }
  else if((r >= R_N) && (r <= R))
    {
      double var = -pow((1.0 - (r/R)),m+1)/r;
      psi_x = var*r_x;
      psi_y = var*r_y;
    }
  else if(r > R)
    {
      psi_x = 0;
      psi_y = 0;   
    }
}

// ********** Divide one cell into two half-sized daughter cells **********
void divide_cell(list<Cell>& Cells_local, std::list<Cell>::iterator &it, Ran& ran, double height){
  double x = (*it).x, y = (*it).y;
  double rN = (*it).N_radius/sqrt(2.0);
  double nx, ny;
  bool set_p = false;
  bool set_n = false;
  bool try_set_cell = true;
  unsigned int n_tries=0;
  Point center;
  center(0) = 0.5*height;
  center(1) = 0.5*height;
  while(try_set_cell){
    n_tries++;
    nx = rN * (2.0*ran.doub() - 1.0);
    ny = sqrt(rN*rN - nx*nx);
    //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    //--------------- Assuming a circular domanin ------
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Point pos_n;
    pos_n(0) = x-nx;
    pos_n(1) = y-ny;
    Point pos_p;
    pos_p(0) = x+nx;
    pos_p(1) = y+ny;
    if(point_distance(pos_p,center)<0.5*height)
      set_p = true;
    else
      set_p = false;
    if(point_distance(pos_n,center)<0.5*height)
      set_n = true;
    else
      set_n = false;
    if((!set_p || !set_n) && n_tries<=20)
      try_set_cell = true;
    else
      try_set_cell = false;
  }
  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //--------------- Place the new cell ---------------
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if(set_p && set_n){
    (*it).set(x+nx, y+ny, rN, (*it).C_radius/sqrt(2.0), (*it).A_radius/sqrt(2.0), (*it).uptake, (*it).time, 5);
    Cell a;  
    a.set(x-nx, y-ny, rN, (*it).C_radius, (*it).A_radius, (*it).uptake, (*it).time, 5);
    Cells_local.insert(it, a);
  }
  else{
    if(set_p)
      (*it).set(x+nx, y+ny, rN, (*it).C_radius/sqrt(2.0), (*it).A_radius/sqrt(2.0), (*it).uptake, (*it).time, 5);
    else if(set_n)
      (*it).set(x-nx, y-ny, rN, (*it).C_radius/sqrt(2.0), (*it).A_radius/sqrt(2.0), (*it).uptake, (*it).time, 5);
    else
      (*it).set(x, y, rN, (*it).C_radius/sqrt(2.0), (*it).A_radius/sqrt(2.0), (*it).uptake, (*it).time, 5);
  }
}

// ********** Set initial conditions in the simulation for number of cells **********
void init_cond_cells(list<Cell>& Cells_local, const Parameters& all_parameters, Ran& ran)
{
  const Real n_radius            = all_parameters.get<Real>("nucleus_radius");
  const Real c_radius            = all_parameters.get<Real>("cell_radius");
  const Real a_radius            = all_parameters.get<Real>("action_radius");
  const Real domain_diameter     = all_parameters.get<Real>("domain_diameter");
  const Real lambda_cell         = all_parameters.get<Real>("lambda_cell");
  const Real initial_con_live    = all_parameters.get<Real>("initial_con_live");
  const Real initial_con_dead    = all_parameters.get<Real>("initial_con_dead");
  const Real tau_P               = all_parameters.get<Real>("cellc_time");
  const Real tau_g1              = all_parameters.get<Real>("g1_time");
  const double nut_ic            = all_parameters.get<Real>("nutrient_ic");
  const double hyp_th            = all_parameters.get<Real>("hypoxic_thrs");
  const double prol_int          = all_parameters.get<Real>("prol_intes");
  //const Real tau_A               = all_parameters.get<Real>("apop_time");
  double confluence_live = 0.0;
  double confluence_dead = 0.0;
  const double proliferative_ratio = max(prol_int*(nut_ic-hyp_th)/(1.0-hyp_th),0.0);
  double single_cell_area = M_PI*std::pow(c_radius,2);
  double domain_area = M_PI*std::pow(0.5*domain_diameter,2);
  // ********** Compute live cells confluence minus 1 cell (to round later) **********
  const double almost_live_confluence = initial_con_live-single_cell_area/domain_area;
  while(confluence_live < almost_live_confluence){
    Point pos;
    bool place_cell = true;
    pos(0) = ran.doub()*domain_diameter;
    pos(1) = ran.doub()*domain_diameter;
    Point center;
    center(0) = 0.5*domain_diameter;
    center(1) = 0.5*domain_diameter;
    if(point_distance(pos, center) > 0.5*domain_diameter)
      place_cell = false;
    else{
      std::list<Cell>::iterator it;
      for(it = Cells_local.begin(); it != Cells_local.end(); it++){
        Point cell;
        cell(0) = (*it).x;
        cell(1) = (*it).y;
        if(point_distance(pos, cell) < 2.*n_radius){
          place_cell = false;
          break;
        }
      }
    }
    if(place_cell){
      Cell a;
      a.set(pos(0), pos(1), n_radius, c_radius, a_radius, lambda_cell, 0, 1);
      if(ran.doub() < proliferative_ratio){
	    a.time = -ran.doub()*(tau_P-tau_g1);
	    a.prev_state = 1;
	    a.state = 2;
	  }
      Cells_local.push_back(a);
      confluence_live += single_cell_area/domain_area;
    }
  }
  // ********** Add missing live cell **********
  double new_radius_cell = sqrt(((initial_con_live-confluence_live)*domain_area)/M_PI);
  double radius_scale = new_radius_cell/c_radius;
  double cell_time = (2.0*pow(radius_scale,2)-1.0)*tau_g1-tau_g1+tau_P;
  bool round_confluence = true;
  while(round_confluence){
    Point pos;
    bool place_cell = true;
    pos(0) = ran.doub()*domain_diameter;
    pos(1) = ran.doub()*domain_diameter;
    Point center;
    center(0) = 0.5*domain_diameter;
    center(1) = 0.5*domain_diameter;
    if(point_distance(pos, center) > 0.5*domain_diameter)
      place_cell = false;
    else{
      std::list<Cell>::iterator it;
      for(it = Cells_local.begin(); it != Cells_local.end(); it++){
        Point cell;
        cell(0) = (*it).x;
        cell(1) = (*it).y;
        if(point_distance(pos, cell) < 2.*n_radius){
          place_cell = false;
          break;
        }
      }
    }
    if(place_cell){
      Cell a;
      a.set(pos(0), pos(1), radius_scale*n_radius, radius_scale*c_radius, radius_scale*a_radius, lambda_cell, 0, 1);
      a.time = -cell_time;
      a.state = 5;
      Cells_local.push_back(a);
      confluence_live += (M_PI*std::pow(radius_scale*c_radius,2))/domain_area;
	  round_confluence = false;
    }
  }
  // ********** Compute dead cells confluence minus 1 cell (to round later) **********
  const double almost_dead_confluence = initial_con_dead-single_cell_area/domain_area;
  //============================== Dead Cells ==============================
  while(confluence_dead < almost_dead_confluence){
    Point pos;
    bool place_cell = true;
    pos(0) = ran.doub()*domain_diameter;
    pos(1) = ran.doub()*domain_diameter;
    Point center;
    center(0) = 0.5*domain_diameter;
    center(1) = 0.5*domain_diameter;
    if(point_distance(pos, center) > 0.5*domain_diameter)
      place_cell = false;
    else{
      std::list<Cell>::iterator it;
      for(it = Cells_local.begin(); it != Cells_local.end(); it++){
        Point cell;
        cell(0) = (*it).x;
        cell(1) = (*it).y;
        if(point_distance(pos, cell) < 2.*n_radius){
          place_cell = false;
          break;
        }
      }
    }
    if(place_cell){
      Cell a;
      a.set(pos(0), pos(1), n_radius, c_radius, a_radius, lambda_cell, 0, 4);
      /*
      a.time = -ran.doub()*(tau_A);
      if(ran.doub() < proliferative_ratio){
	    a.time = -ran.doub()*(tau_P-tau_g1);
	    a.prev_state = 1;
	    a.state = 2;
	  }
	  */
      Cells_local.push_back(a);
      confluence_dead += ((M_PI*std::pow(c_radius,2))/(M_PI*std::pow(0.5*domain_diameter,2)));
    }
  }
  // ********** Add missing dead cell **********
  new_radius_cell = sqrt(((initial_con_dead-confluence_dead)*domain_area)/M_PI);
  radius_scale = new_radius_cell/c_radius;
  cell_time = (2.0*pow(radius_scale,2)-1.0)*tau_g1-tau_g1+tau_P;  
  round_confluence = true;
  while(round_confluence){
    Point pos;
    bool place_cell = true;
    pos(0) = ran.doub()*domain_diameter;
    pos(1) = ran.doub()*domain_diameter;
    Point center;
    center(0) = 0.5*domain_diameter;
    center(1) = 0.5*domain_diameter;
    if(point_distance(pos, center) > 0.5*domain_diameter)
      place_cell = false;
    else{
      std::list<Cell>::iterator it;
      for(it = Cells_local.begin(); it != Cells_local.end(); it++){
        Point cell;
        cell(0) = (*it).x;
        cell(1) = (*it).y;
        if(point_distance(pos, cell) < 2.*n_radius){
          place_cell = false;
          break;
        }
      }
    }
    if(place_cell){
      Cell a;
      a.set(pos(0), pos(1), radius_scale*n_radius, radius_scale*c_radius, radius_scale*a_radius, lambda_cell, 0, 4);
      /*
	  a.time = -ran.doub()*(tau_A);
      if(ran.doub() < proliferative_ratio){
	    a.time = -ran.doub()*(tau_P-tau_g1);
	    a.prev_state = 1;
	    a.state = 2;
	  }
	  */
      Cells_local.push_back(a);
      confluence_dead += (M_PI*std::pow(radius_scale*c_radius,2))/domain_area;
	  round_confluence = false;
    }
  }
}

// ********** Calculate distance between a cell and a point **********
double distance(Cell cell_A, double pos_x, double pos_y)
{
  return sqrt(pow(cell_A.x-pos_x,2)+pow(cell_A.y-pos_y,2));
}

// ********** Save cells for output **********
void save_cells(const list<Cell>& Cells_local, double domain_diameter, string s, int file_number, int t)
{
  bool matlab = true;
  if(matlab){
    //==******** Creating a string with the correct name ********==//
    const char *c = s.c_str();
    char n[100],name[200];
    sprintf(n,"%d_%05d.m",file_number,t);
    strcpy(name,c);
    strcat(name,n);
    stringstream ss;
    string name_s;
    ss << name;
    ss >> name_s;

    //==******** Saving the data ********==//
    ofstream out_file;
    out_file.open (name_s);
    out_file << "cells = zeros(" << Cells_local.size()+1 << "," << 4 << ");" <<endl;
    out_file << "cells = [" << -1 << " " << domain_diameter*0.5 << " " << domain_diameter*0.5 << " " << domain_diameter*0.5 << endl;
  
    // Save all cells
    std::list<Cell>::const_iterator it;
    for(it = Cells_local.begin(); it != Cells_local.end(); it++)
    {
      out_file << (*it).state << " " << (*it).x << " " << (*it).y << " " << (*it).C_radius << endl;
    }
  
    out_file << "];";
    out_file.close();
  }
  else{
    //========== Creating a string with the correct name ==========//
    std::list<Cell>::const_iterator it;
    unsigned int dead = 0;
    unsigned int live = 0;
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      if ( (*it).state==4 )
        dead++;
      else
        live++;
    }
    const char *c = s.c_str();
    char n[100],name[200];
    sprintf(n,"%d_%05d.txt",file_number,t);
    strcpy(name,c);
    strcat(name,n);
    stringstream ss;
    string name_s;
    ss << name;
    ss >> name_s;
    //========== Saving the data ==========//
    ofstream out_file;
    out_file.open (name_s);
    out_file << domain_diameter << " " << domain_diameter << endl;
    out_file << Cells_local.size() << " " << t << endl;
    out_file << live << " " << dead << endl;
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      out_file << (*it).state << endl;
      out_file << scientific;
      out_file << (*it).x << " " << (*it).y << endl;
      out_file << (*it).N_radius << " " << (*it).C_radius << " " << (*it).cal << endl;
      out_file << (*it).v[0] << " " << (*it).v[1] << endl;
    }  
    out_file.close();
  }
}

// ********** Cell member function to print state values **********
void Cell::print()
{
  cout << "Position = ( " << x << " , " << y << " )" << endl;
  cout << "Radius | Nuclear = " << N_radius << ", Cell = " << C_radius << ", Max = " << A_radius << endl;
  cout << "Time = " << time << endl;
  cout << "V = (" << v[0] << " , " << v[1] << " )" << endl;
  
  switch (state) 
    {
    case 0:
      cout << "State = Dead" << endl;
      break;
    case 1:
      cout << "State = Alive" << endl;
      break;
    case 2:
      cout << "State = Proliferative" << endl;
      break;
    case 3:
      cout << "State = Hypoxic" << endl;
      break;
    case 4:
      cout << "State = Dying" << endl;
      break;
    case 5:
      cout << "State = G1" << endl;
      break;
    case 6:
      cout << "State = Normoxic" << endl;
      break;
    default:
      cout << "State = " << state << endl;
    }
}

// ********** Set member values of cell object **********
void Cell::set(double X,double Y,double RN,double R,double RA,double Uptake, int Time,int State)
{
  x = X;
  y = Y;
  N_radius = RN;
  C_radius = R;
  A_radius = RA;
  uptake = Uptake;
  time = Time;
  state = State;
  v[0] = 0.0;
  v[1] = 0.0;
  F[0] = 0.0;
  F[1] = 0.0;
  cal  = 0.0;
}
