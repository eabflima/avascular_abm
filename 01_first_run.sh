#!/bin/bash
rm *# *~ &> /dev/null
scale=16
#################### Create Options File ######################
rm *# *~ options.in &> /dev/null
echo "mesh_name = circle_6384_lc400.msh" >> options.in
echo "read_mesh = 1"                     >> options.in
echo "ic_type = 27"                      >> options.in
echo "initial_tum = 3"                   >> options.in
echo "n_timesteps = 96"                  >> options.in
echo "print_inter = 3"                   >> options.in
echo "verbose = 0"                       >> options.in
echo "print_sa = 0"                      >> options.in
echo "domain_diameter = 6384"            >> options.in
echo "cluster_scale = ${scale}"          >> options.in
echo "con_b = 0.001"                     >> options.in
echo "con_n = 0.001"                     >> options.in
echo "time_step = 1.0"                   >> options.in
echo "c_ccr = 10.0"                      >> options.in
echo "c_tta = 0.488836"                  >> options.in
echo "c_hha = 0.588836"                  >> options.in
echo "max_outside = 2000"                >> options.in
echo "rand_seed = 5"                     >> options.in
echo "prol_intens = 0.52"                >> options.in
echo "file_number = 1"                   >> options.in
echo "k_var = 1"                         >> options.in
echo "k_val = 0.1"                       >> options.in
echo "nucleus_radius = 5.295"            >> options.in
echo "cell_radius = 9.953"               >> options.in
echo "action_prop = 1.214"               >> options.in
echo "ntri_ic = 0.5"                     >> options.in
#################### Create Python Mean Var Code ##############
rm *# *~ mean_std.py &> /dev/null
echo "import sys"                                              >> mean_std.py
echo "import pandas as pd"                                     >> mean_std.py
echo "data = pd.read_csv('sbl.csv',header=None,delimiter=' ')" >> mean_std.py
echo "orig_stdout = sys.stdout"                                >> mean_std.py
echo "f = open('out_py.txt', 'w')"                             >> mean_std.py
echo "sys.stdout = f"                                          >> mean_std.py
echo "for column in data:"                                     >> mean_std.py
echo "    print data[column].mean(),data[column].std()"        >> mean_std.py
echo "sys.stdout = orig_stdout"                                >> mean_std.py
echo "f.close()"                                               >> mean_std.py
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Run the model --------------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
loop=40
echo "Realizations to add = ${loop}"
#################### Create Parameters File ###################
rm *# *~ parameters.in &> /dev/null
echo "alpha_p = 4.9e-02"   >> parameters.in
echo "alpha_a = 4.1e-04"   >> parameters.in
echo "nu_diff = 50.0"      >> parameters.in
echo "n_lambd = 4.8e-02"   >> parameters.in
echo "live_ic = 0.5"       >> parameters.in
echo "dead_ic = 0.3"       >> parameters.in
echo "rate_pc = 0.0"       >> parameters.in
echo "t_death = 97.0"      >> parameters.in
echo "gamma_a = 2.4e-02"   >> parameters.in
echo "gamma_p = 50.0"      >> parameters.in
echo "sigma_h = 5.4e-02"   >> parameters.in
echo "size_loop = ${loop}" >> parameters.in
FILE=sbl.csv
if [ -f "$FILE" ]; then
    cp $FILE tmp.csv
    make run
    cat tmp.csv >> sbl.csv
else 
    make run
fi
number=$(more $FILE | wc -l | awk '{printf "%05d\n",$1+0}')
echo $number
pvpython mean_std.py
sed -n -e "1,33p" out_py.txt > live_${number}.dat
sed -n -e "34,66p" out_py.txt > dead_${number}.dat
rm *# *~ tmp.csv out_py.txt &> /dev/null
