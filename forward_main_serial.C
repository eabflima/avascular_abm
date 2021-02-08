// Running the main function as a function to test for vegf_prod values.
#include "general_libraries.h"
#include "main_abm_oxy.h"
#include "cell_abm.h"
#include "mpi.h" 

int main(int argc, char** argv){
  
  MPI_Init(&argc, &argv);
  
  {
    LibMeshInit init (argc, argv, MPI_COMM_SELF);
    
    std::vector<double> Parameters(11, 0);
    
    GetPot input_file("parameters.in");
    Parameters[0]  = input_file("alpha_p",1.0);
    Parameters[1]  = input_file("alpha_a",1.0);
    Parameters[2]  = input_file("nu_diff",1.0);
    Parameters[3]  = input_file("n_lambd",1.0);
    Parameters[4]  = input_file("live_ic",1.0);
    Parameters[5]  = input_file("dead_ic",1.0);
    Parameters[6]  = input_file("rate_pc",1.0);
    Parameters[7]  = input_file("t_death",1.0);
    Parameters[8]  = input_file("gamma_a",1.0);
    Parameters[9]  = input_file("gamma_p",1.0);
    Parameters[10] = input_file("sigma_h",1.0);
       
    const unsigned int size = input_file("size_loop",1);
    srand(time(NULL));
    std::vector<unsigned int> seeds(size,0);
    for(unsigned int i = 0; i < size; i++)
      seeds[i] = rand();
    ofstream out_file;
    out_file.precision(16);
    out_file.open("sbl.csv");
    for(unsigned int j=0;j<size;j++){
      cout << "Seed number = " << j << endl;
      std::vector<double> liveConfModel;
      std::vector<double> deadConfModel;
      main_code(init, liveConfModel, deadConfModel, Parameters, seeds[j]);
      out_file << scientific << liveConfModel[0];
      for(unsigned int o=1; o<liveConfModel.size(); o++)
	out_file << " " << liveConfModel[o];
      for(unsigned int o=0; o<deadConfModel.size(); o++)
	out_file << " " << deadConfModel[o];
      out_file << endl;
    }
    out_file.close();
  }
  MPI_Finalize();
}
