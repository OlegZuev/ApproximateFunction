#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>

#include "Functions.h"
constexpr int COUNT_OF_ROWS = N + 1;
constexpr int COUNT_OF_COLUMNS = N + 2;

using namespace std;

const bool PRINT_TO_FILE = false;



int main() {
	ofstream fout("../" + to_string(VARIANT) + "_output" + ".txt");
	ostream& out = PRINT_TO_FILE ? fout : cout;
	
	double** table_SP;
	allocate_matrix(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS);
	for (int i = 0; i < COUNT_OF_ROWS; ++i) {
		table_SP[i][0] = X_0 + i * (X_N - X_0) / 5;
	}

	calculate_table_of_split_differences(table_SP, COUNT_OF_ROWS, COUNT_OF_ROWS, func);
	print_matrix(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS, out);

	newton_method(table_SP, out);
	spline_interpolation_method(table_SP, out);
	mean_square_approximation_discrete(table_SP, out);
	mean_square_approximation_integral(table_SP, out);


	delete_matrix(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS);
}
