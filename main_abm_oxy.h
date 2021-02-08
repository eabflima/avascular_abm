#ifndef MAIN_ABM_OXY_H_
#define MAIN_ABM_OXY_H_

#include "general_libraries.h"
#include "cell_abm.h"

Number exact_value_ox(const Point&, const Parameters&, const std::string&, const std::string& );

double triangle_area(const Point& e1, const Point& e2, const Point& e3);

double scalar_prod(const Point& e1, const Point& e2);

void initial_condition_ox(EquationSystems& es, const std::string& system_name);

void assemble_ox(EquationSystems& es,  const std::string& system_name);

double area_inside_tri(const Point& p, const Point& e1, const Point& e2, const Point& e3, double radius, double elem_volume);

double area_inside_quad(const Point& p, const Point& e1, const Point& e2, const Point& e3, const Point& e4, double radius);

void main_code(LibMeshInit &init, vector<double>& LIVE_CONFL, vector<double>& DEAD_CONFL, vector<double> Parameters, const unsigned int rand_seed);

#endif
