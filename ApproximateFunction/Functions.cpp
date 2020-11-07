#include "Functions.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

/*
 * table_SP - table_of_split_differences.
 * n - count of rows of actual matrix
 * k - count of computing elements in column
 */
void calculate_table_of_split_differences(double** table_SP, int n, int k, double (*func)(double)) {
	if (n - k == 0) {
		for (int i = 0; i < k; ++i) {
			table_SP[i][n - k + 1] = func(table_SP[i][n - k]);
		}
	} else {
		for (int i = 0; i < k; ++i) {
			table_SP[i][n - k + 1] = (table_SP[i + 1][n - k] - table_SP[i][n - k]) / (table_SP[n - k][0] - table_SP[0][0
			]);
		}
	}

	if (k - 1 != 0) {
		calculate_table_of_split_differences(table_SP, n, k - 1, func);
	}
}

/*
 * func for task
 */
double func(double x) {
	switch (VARIANT) {
	case 0:
		return pow(3, x - 1) + 4 - x;
	case 5:
		return pow(3, x - 1) + 2 - x;
	case 7:
		return exp(-2 * x) - 2 * x + 1;
	}

	throw std::runtime_error("Incorrect variant");
}

/*
 * approximate func computed from root mean square approximation
 */
double func_mean_square(double x, double* c) {
	return c[0] + c[1] * x + c[2] * pow(x, 2);
}

/*
 * max f^(number_of_derivative) on [a,b]
 */
double compute_m(int number_of_derivative) {
	switch (VARIANT) {
	case 0:
	case 5:
		return 3 * pow(log(3), number_of_derivative);
	case 7:
		return pow(-2, number_of_derivative) * exp(-4);
	}

	throw std::runtime_error("I don't know how to compute this derivative");
}

/*
 *
 */
double compute_derivative_func(double x) {
	switch (VARIANT) {
	case 0:
	case 5:
		return pow(3, x - 1) * log(3) - 1;
	case 7:
		return -2 * exp(-4) - 2;
	}

	throw std::runtime_error("Incorrect variant");
}

double func_Pn(double** table_SP, double x, double q, int k) {
	if (k >= N + 1) {
		return 0;
	}

	if (k == 0) {
		q = 1;
	} else {
		q *= x - table_SP[k - 1][0];
	}

	return table_SP[0][k + 1] * q + func_Pn(table_SP, x, q, k + 1);
}

double find_error_for_newton(double** table_SP, double x, int n, double M6) {
	double omega = 1;
	for (int i = 0; i < n + 1; i++) {
		omega *= x - table_SP[i][0];
	}

	int fact = 1;
	for (int i = 2; i < n + 2; i++) {
		fact *= i;
	}

	return M6 * fabs(omega) / fact;
}

double find_error_for_spline_interpolation_method(double M4, double M5) {
	return (M4 / 384 + M5 * H / 240) * pow(H, 4);
}

void newton_method(double** table_SP, std::ostream& ostr) {
	double M6 = compute_m(6);
	ostr << "M6 = " << M6 << std::endl;
	for (int i = 0; i < N; ++i) {
		double x = X_0 + (i + 0.5) * (X_N - X_0) / 5;
		ostr << x << " " << func(x) << " " << func_Pn(table_SP, x, 1, 0) << " " <<
			fabs(func(x) - func_Pn(table_SP, x, 1, 0)) << " " << find_error_for_newton(table_SP, x, N, M6) << std::endl;
	}
}

void spline_interpolation_method(double** table_SP, std::ostream& ostr) {
	double m[N + 1];
	m[0] = compute_derivative_func(X_0);
	m[N] = compute_derivative_func(X_N);
	double alpha[N + 1];
	double beta[N + 1];
	alpha[1] = 0;
	beta[1] = m[0];
	// alpha[1] = -0.25;
	// beta[1] = 3 * (func(table_SP[1][0]) - func(0.8)) / H / 4;
	//beta[1] = (3 * (func(table_SP[2][0]) - func(table_SP[0][0])) / H - m[0]) / 4;

	for (int j = 1; j < N; j++) {
		alpha[j + 1] = -1 / (4 + alpha[j]);
		beta[j + 1] = (3 * (func(table_SP[j + 1][0]) - func(table_SP[j - 1][0])) / H - beta[j]) / (4 + alpha[j]);
	}

	for (int j = N - 1; j >= 0; j--) {
		m[j] = alpha[j + 1] * m[j + 1] + beta[j + 1];
	}

	double M5 = compute_m(5);
	ostr << "M5 = " << M5 << std::endl << std::endl;

	ostr << " x[i] | " << " df/dx(x[i]) " << "  m[i]   " << " Delta " << " Estimate " << std::endl;
	for (int i = 0; i < N + 1; i++) {
		ostr << std::setprecision(2) << table_SP[i][0] << "   | " << std::setprecision(6) <<
			compute_derivative_func(table_SP[i][0]) << " | " << m[i] << " | " << fabs(
				compute_derivative_func(table_SP[i][0]) - m[i]) << " | " << M5 / 60 * pow(H, 4) << std::endl;
	}

	double M4 = compute_m(4);
	ostr << "M4 = " << M4 << std::endl << std::endl;

	double x_arr[N + 1];
	for (int i = 0; i < N + 1; i++) {
		x_arr[i] = X_0 + (i + 0.5) * (X_N - X_0) / 5;
	}

	ostr << std::endl << " x  | " << "   f(x)   " << "   S(f, x) " << " Delta " << " Estimate " << std::endl;
	for (int i = 0; i < N; i++) {
		double tau = (x_arr[i] - table_SP[i][0]) / H;
		double S = fi0(tau) * func(table_SP[i][0]) + fi0(1 - tau) * func(table_SP[i + 1][0]) + H * (fi1(tau) * m[i] -
			fi1(1 - tau) *
			m[i + 1]);
		ostr << std::setprecision(2) << x_arr[i] << "  | " << std::setprecision(6) << func(x_arr[i]) <<
			" | " << S << " | " << abs(S - func(x_arr[i])) << " | " <<
			find_error_for_spline_interpolation_method(M4, M5) << std::endl;
	}
}

void mean_square_approximation_discrete(double** table_SP, std::ostream& ostr) {
	ostr << "Discrete variant" << std::endl << std::endl;

	double** matrix;
	allocate_matrix(matrix, 3, 3);

	matrix[0][0] = N + 1;
	for (int i = 0; i < N + 1; i++) {
		matrix[1][0] += table_SP[i][0];
		matrix[2][0] += pow(table_SP[i][0], 2);
		matrix[2][1] += pow(table_SP[i][0], 3);
		matrix[2][2] += pow(table_SP[i][0], 4);
	}

	matrix[0][1] = matrix[1][0];
	matrix[0][2] = matrix[2][0];
	matrix[1][1] = matrix[2][0];
	matrix[1][2] = matrix[2][1];

	ostr << "Matrix:" << std::endl;
	print_matrix(matrix, 3, 3, ostr);

	double* b = new double[3]{0};
	for (int i = 0; i < N + 1; ++i) {
		b[0] += func(table_SP[i][0]);
		b[1] += func(table_SP[i][0]) * table_SP[i][0];
		b[2] += func(table_SP[i][0]) * pow(table_SP[i][0], 2);
	}

	ostr << "Vector of right parts:" << std::endl;
	for (int i = 0; i < 3; ++i) {
		ostr << b[i] << " ";
	}

	ostr << std::endl;
	double* c = gauss_method(matrix, b, 3);

	ostr << std::setprecision(4) << "P2(x) = " << c[0] << " + (" << c[1] << "*x) + (" << c[2] << "*x^2)\n";
	double sum = 0;
	for (int i = 0; i < N + 1; i++) {
		sum += pow(func(table_SP[i][0]), 2) - pow(func_mean_square(table_SP[i][0], c), 2);
	}

	double error = sqrt(sum);
	ostr << "Error estimation: " << error << std::endl;
	delete_matrix(matrix, 3, 3);
	delete[] b;
	delete[] c;
}

void mean_square_approximation_integral(double** table_SP, std::ostream& ostr) {
	ostr << "Integral variant" << std::endl << std::endl;

	double** matrix;
	allocate_matrix(matrix, 3, 3);
	// integrals of the form (g, g) = ∫g(x)g(x)dx, x = 1..2 where g1 = 1, g2 = x, g3 = x^2
	matrix[0][0] = 1.0;
	matrix[0][1] = 3.0 / 2;
	matrix[0][2] = 7.0 / 3;
	matrix[1][0] = 3.0 / 2;
	matrix[1][1] = 7.0 / 3;
	matrix[1][2] = 15.0 / 4;
	matrix[2][0] = 7.0 / 3;
	matrix[2][1] = 15.0 / 4;
	matrix[2][2] = 31.0 / 5;

	ostr << "Matrix:" << std::endl;
	print_matrix(matrix, 3, 3, ostr);

	double* b = new double[3]{0};
	// integrals of the form (f, g) = ∫f(x)g(x)dx, x = 1..2 where g1 = 1, g2 = x, g3 = x^2
	switch (VARIANT) {
	case 0:
		b[0] = 4.32047845325358;
		b[1] = 6.56079190041974;
		b[2] = 10.3272721971579;
		break;
	case 5:
		b[0] = (4 + log(3)) / log(9);
		b[1] = 2.0 / 3 + (log(243) - 2) / pow(log(3), 2);
		b[2] = 11.0 / 12 + (4 + 11 * pow(log(3), 2) - 5 * log(9)) / pow(log(3), 3);
		break;
	case 7:
		b[0] = -2 + (-1 + exp(2)) / (2 * exp(4));
		b[1] = -19.0 / 6 + (-5 + 3 * exp(2)) / (4 * exp(4));
		b[2] = -31.0 / 6 + (-13 + 5 * exp(2)) / (4 * exp(4));
		break;
	default:
		throw std::runtime_error("Incorrect variant");
	}

	ostr << "Vector of right parts:" << std::endl;
	for (int i = 0; i < 3; ++i) {
		ostr << b[i] << " ";
	}

	ostr << std::endl;
	double* c = gauss_method(matrix, b, 3);

	ostr << std::setprecision(4) << "P2(x) = " << c[0] << " + (" << c[1] << "*x) + (" << c[2] << "*x^2)\n";

	double ff;
	double fg;
	switch (VARIANT) {
	case 0:
		ff = 0.0000543668793042684705375205705609;
		fg = 0;
		break;
	case 5:
		ff = 1.0 / 3 + (4 + log(9)) / pow(log(3), 2);
		fg = (c[2] * (4 + 11 * pow(log(3), 2) - 5 * log(9)) + log(3) * (c[0] * log(9) + c[1] * (-2 + log(243)))) /
			pow(log(3), 3) + 1.0 / 12 * (6 * c[0] + 8 * c[1] + 11 * c[2]);
		break;
	case 7:
		ff = 13.0 / 3 - 1 / (4 * exp(8)) + 17 / (4 * exp(4)) - 2 / exp(2);
		fg = (c[2] * (-39 + 15 * exp(2) - 62 * exp(4)) + c[1] * (-15 + 9 * exp(2) - 38 * exp(4)) - 6 * c[0] * (1 -
			exp(2) + 4 * exp(4))) / (12 * exp(4));
		break;
	default:
		throw std::runtime_error("Incorrect variant");
	}

	double error = sqrt(ff - fg);
	ostr << "Error estimation: " << error << std::endl;
	delete_matrix(matrix, 3, 3);
	delete[] b;
	delete[] c;
}

double* gauss_method(double** matrix, double* y, int n) {
	double* x = new double[n]{0};
	const double eps = 0.00001;
	for (int k = 0; k < n; k++) {
		// find max elem in column
		double max = fabs(matrix[k][k]);
		int index = k;
		for (int i = k + 1; i < n; i++) {
			if (fabs(matrix[i][k]) > max) {
				max = fabs(matrix[i][k]);
				index = i;
			}
		}

		if (max < eps) {
			std::cout << "No solution";
			return nullptr;
		}

		// Permutation of rows
		std::swap(matrix[k], matrix[index]);
		std::swap(y[k], y[index]);

		// Normalization of equations
		for (int i = k; i < n; i++) {
			double multiplier = matrix[i][k];
			if (fabs(multiplier) < eps) {
				continue;
			}

			for (int j = 0; j < n; j++) {
				matrix[i][j] = matrix[i][j] / multiplier;
			}

			y[i] = y[i] / multiplier;
			if (i == k) {
				continue;
			}

			for (int j = 0; j < n; j++) {
				matrix[i][j] = matrix[i][j] - matrix[k][j];
			}

			y[i] = y[i] - y[k];
		}
	}

	// reverse substitution
	for (int k = n - 1; k >= 0; k--) {
		x[k] = y[k];
		for (int i = 0; i < k; i++) {
			y[i] = y[i] - matrix[i][k] * x[k];
		}
	}

	return x;
}

/*
 * allocate matrix n * m
 * n - count of rows
 * m - count of columns
 */
void allocate_matrix(double**& matrix, int n, int m) {
	matrix = new double*[n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new double[m]{0};
	}
}

/*
 * delete matrix n * m
 * n - count of rows
 * m - count of columns
 */
void delete_matrix(double** matrix, int n, int m) {
	for (int i = 0; i < n; ++i) {
		delete[] matrix[i];
	}

	delete[] matrix;
}

void print_matrix(double** matrix, int n, int m, std::ostream& ostr) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			ostr << matrix[i][j] << " ";
		}
		ostr << std::endl;
	}
	ostr << std::endl;
}

double fi0(double tau) {
	return (1 + 2 * tau) * pow(1 - tau, 2);
}

double fi1(double tau) {
	return tau * pow(1 - tau, 2);
}
