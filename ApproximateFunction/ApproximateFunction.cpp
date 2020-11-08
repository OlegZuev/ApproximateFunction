#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include "Functions.h"
using namespace std;

const bool PRINT_TO_FILE = true;

int main() {
	ofstream fout("../" + to_string(VARIANT) + "_output" + ".txt");
	ostream& out = PRINT_TO_FILE ? fout : cout;
	
	double** table_SP;
	allocate_matrix(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS);
	for (int i = 0; i < COUNT_OF_ROWS; ++i) {
		table_SP[i][0] = X_0 + i * (X_N - X_0) / 5;
	}

	out << "Variant: " << VARIANT << endl;
	calculate_table_of_split_differences(table_SP, COUNT_OF_ROWS, COUNT_OF_ROWS, func);
	out << "Newton's interpolation formula" << endl;
	print_table_of_split_differences(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS, out, 2);

	newton_method(table_SP, out);
	spline_interpolation_method(table_SP, out);
	out << "Mean square approximation" << endl << endl;
	mean_square_approximation_discrete(table_SP, out);
	mean_square_approximation_integral(table_SP, out);
	reverse_interpolation_method(table_SP, out);

	delete_matrix(table_SP, COUNT_OF_ROWS, COUNT_OF_COLUMNS);
}
