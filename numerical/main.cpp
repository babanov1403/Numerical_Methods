#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "integrate_methods.h"
#include "matrix_class.h"
#include "burger.h"
#include "polynom_class.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"


using namespace biv;
using namespace std;
int main() {
	Matrix<double> RawC = buildC(t0, t1, N, A, B, H, g, x0);
	Matrix<double> RawD = buildD(t0, t1, N, A, B, H, g, x0);
	Matrix<double> RawG0 = compute_g0(g, H, x0);

	vector<double> U(N);
	alglib::minlpstate state;
	alglib::minlpreport rep;
	alglib::real_1d_array x;
	solveForCD(state, rep, x, RawC, RawD, RawG0);
	cout << "x:\n";
	
	for (int i = 0; i < N; i++) U[i] = x(i);

	for (int i = 0; i < N; i++) cout << U[i] << ' ';

	cout << "\nTarget:\n";
	double ans = 0;
	for (int i = 0; i < N; i++) 
		ans += RawC(0, i) * x(i);
	cout << -ans << '\n';

	cout << "Computing:\n";
	Matrix<double> FinalPos = computePathWGivenU(U, true);
	cout << "FinalPos:\n";
	cout << FinalPos;

	Matrix<double> DoesIt = H * FinalPos;
	cout << "Nevyazka:\n";

	cout << DoesIt << g;

}
