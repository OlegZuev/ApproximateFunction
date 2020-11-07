#pragma once
#include <ostream>

const int VARIANT = 0;
const int N = 5;
const double X_0 = 1;
const double X_N = 2;
const double H = 0.2;

void calculate_table_of_split_differences(double** table_SP, int n, int k, double (*func)(double));

void spline_interpolation_method(double** table_SP, std::ostream& ostr);

void mean_square_approximation_discrete(double** table_SP, std::ostream& ostr);

void mean_square_approximation_integral(double** table_SP, std::ostream& ostr);

double func(double x);

double compute_m(int number_of_derivative);

double find_error_for_newton(double** table_SP, double x, int n, double M6);

void newton_method(double** table_SP, std::ostream& ostr);

void allocate_matrix(double**& matrix, int n, int m);

void delete_matrix(double** matrix, int n, int m);

void print_matrix(double** matrix, int n, int m, std::ostream& ostr);

double func_Pn(double** table_SP, double x, double q, int k);

double fi0(double tau);

double fi1(double tau);

double* gauss_method(double** matrix, double* y, int n);