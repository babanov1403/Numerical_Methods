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
	Matrix<double> RawC = buildC(t0, t1, N, A, B, H, g, x0);
	cout << "C:\n";
	cout << RawC;
	Matrix<double> RawD = buildD(t0, t1, N, A, B, H, g, x0);
	cout << "D:\n";
	cout << RawD;
	cout << "g0:\n";
	Matrix<double> RawG0  = compute_g0(g, H, x0);
	cout << RawG0;
	//у нас есть вектор C 1xN, вектор D mxN, g0 mx1, и граниченые условия |u_i| < 1 =>>

	alglib::real_2d_array D; D.setlength(H.n, N);
	alglib::real_1d_array Target; Target.setlength(N);
	alglib::real_1d_array LBoundaries_u; LBoundaries_u.setlength(N);
	alglib::real_1d_array UBoundaries_u; UBoundaries_u.setlength(N);
	alglib::real_1d_array LBoundaries_D; LBoundaries_D.setlength(H.n);
	alglib::real_1d_array UBoundaries_D; UBoundaries_D.setlength(H.n);
	alglib::real_1d_array Scale; Scale.setlength(N);
	alglib::minlpstate state;
	alglib::minlpreport rep;
	alglib::real_1d_array x;
	alglib::minlpcreate(N, state);
	//filling arrays
	for (int i = 0; i < RawD.n; i++) 
		for (int j = 0; j < RawD.m; j++)
			D(i, j) = RawD(i, j);
	
	for (int i = 0; i < RawD.n; i++){
		for (int j = 0; j < RawD.m; j++)
			cout << D(i, j) << ' ';
		cout << "\n";
	}

	for (int i = 0; i < N; i++)
		Target(i) = -RawC(0, i);
	
	for (int i = 0; i < N; i++)
		cout << Target(i) << ' ';
	cout << '\n';

	for (int i = 0; i < N; i++) {
		LBoundaries_u(i) = -10;
		UBoundaries_u(i) = 10;
	} 

	for (int i = 0; i < N; i++) {
		cout << LBoundaries_u(i)  << ' ';
		cout << UBoundaries_u(i)  << ' ';
		cout << '\n';
	}
	
	for (int i = 0; i < H.n; i++) {
		LBoundaries_D(i) = RawG0(i, 0);
		UBoundaries_D(i) = LBoundaries_D(i);
	}
	for (int i = 0; i < H.n; i++) {
		cout << LBoundaries_D(i) << ' ';
		cout << UBoundaries_D(i) << ' ';
		cout << '\n';
	}

	for (int i = 0; i < N; i++)
		Scale(i) = 1;
	//////
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
	for (int i = 0; i < N; i++) cout << x(i) << ' ';
	cout << "\nTarget:\n";
	double ans = 0;
	for (int i = 0; i < N; i++) 
		ans += Target(i) * x(i);
	cout << -ans << '\n';

	vector<double> U(N);
	for (int i = 0; i < N; i++) U[i] = x(i);
	cout << "Computing:\n";
	Matrix<double> FinalPos = computePathWGivenU(U);
	cout << "FinalPos:\n";
	cout << FinalPos;
	Matrix<double> DoesIt = H * FinalPos;
	cout << "Nevyazka:\n";
	cout << DoesIt << g;

}
