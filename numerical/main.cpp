#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "integrate_methods.h"
#include "matrix_class.h"
#include "burger.h"
#include "polynom_class.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"


using namespace biv;
using namespace std;
int main() {
	alglib::real_2d_array D = "[[-2, -1, 2, 0, 0], [1, 1, 4, 1, 3], [3, 1, -1, 0, 6]]";
	alglib::real_1d_array Target = "[0, 1, -3, -1, -1]";
	alglib::real_1d_array LBoundaries_u = "[0, 0, 0, 0, 0]";
	alglib::real_1d_array UBoundaries_u = "[+inf, +inf, +inf, +inf, +inf]";
	alglib::real_1d_array LBoundaries_D = "[2, 8, 5]";
	alglib::real_1d_array UBoundaries_D = "[2, 8, 5]";
	alglib::real_1d_array Scale = "[1, 1, 1, 1, 1]";
	alglib::minlpstate state;
	alglib::minlpreport rep;
	alglib::real_1d_array x;
	alglib::minlpcreate(5, state);
	
	alglib::minlpsetcost(state, Target);
	alglib::minlpsetbc(state, LBoundaries_u, UBoundaries_u);
	alglib::minlpsetlc2dense(state, D, LBoundaries_D, UBoundaries_D);
	alglib::minlpsetscale(state, Scale);

	//solving
	alglib::minlpoptimize(state);
	alglib::minlpresults(state, x, rep);
	cout << "State:\n";
	cout << rep.terminationtype << '\n';
	cout << "x:\n";
	for (int i = 0; i < 5; i++) cout << x(i) << ' ';

}
