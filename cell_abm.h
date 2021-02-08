#ifndef CELL_ABM_H_
#define CELL_ABM_H_

class Cell{
public:
  //********** Cell position **********
  double x, y;
  //********** Cell properties **********
  double N_radius, C_radius, A_radius, uptake, cal;
  double v[2], F[2];
  int time, state, prev_state;
  void set(double X, double Y, double RN, double R, double RA, double Uptake, int Time, int State);
  void print();
};

void init_cond_cells(list<Cell>& Cells_local, const Parameters& all_parameters, Ran& ran);

void save_cells(const list<Cell>& Cells_local,double domain_diameter, string s, int file_number, int t);

double distance(Cell cell_A, double pos_x, double pos_y);

void divide_cell(list<Cell>& Cells_local, std::list<Cell>::iterator &it, Ran& ran, double height);

double euclid_norm(double x, double y);

void potential_adh(double& phi_x, double& phi_y, double r_x, double r_y, double R_A, int n);

void potential_rep(double& psi_x, double& psi_y, double r_x, double r_y, double R_N, double R, double M, int m);

void normal(double& N_x, double& N_y, double x_cel, double y_cel, double height);

void compute_forces(list<Cell>& Cells_local, EquationSystems& equation_systems, double height, int& outside_cells, int& total_tumor);

void update_states(list<Cell> & Cells_local, int time, Ran & ran, EquationSystems& equation_systems, const int outside_cells, int & total_tumor, int & num_dead);

double point_distance(const Point& A, const Point& B);

#endif
