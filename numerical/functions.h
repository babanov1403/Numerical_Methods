#pragma once
#include "optimization.h"
#include "burger.h"
using namespace biv;
namespace biv {

	template<typename T>
	void solveForCD(alglib::minlpstate& state
		, alglib::minlpreport& rep
		, alglib::real_1d_array& x
		, const Matrix<T>& C
		, const Matrix<T>& D
		, const Matrix<T>& G0) {
		alglib::real_2d_array D1; D1.setlength(H.n, N);
		alglib::real_1d_array Target; Target.setlength(N);
		alglib::real_1d_array LBoundaries_u; LBoundaries_u.setlength(N);
		alglib::real_1d_array UBoundaries_u; UBoundaries_u.setlength(N);
		alglib::real_1d_array LBoundaries_D; LBoundaries_D.setlength(H.n);
		alglib::real_1d_array UBoundaries_D; UBoundaries_D.setlength(H.n);
		alglib::real_1d_array Scale; Scale.setlength(N);

		alglib::minlpcreate(N, state);
		for (int i = 0; i <D.n; i++)
			for (int j = 0; j < D.m; j++)
				D1(i, j) = D(i, j);

		for (int i = 0; i < N; i++)
			Target(i) = C(0, i);

		for (int i = 0; i < N; i++) {
			LBoundaries_u(i) = -1;
			UBoundaries_u(i) = 1;
		}

		for (int i = 0; i < H.n; i++) {
			LBoundaries_D(i) = G0(i, 0);
			UBoundaries_D(i) = LBoundaries_D(i);
		}

		for (int i = 0; i < N; i++)
			Scale(i) = 1;

		alglib::minlpsetcost(state, Target);
		alglib::minlpsetbc(state, LBoundaries_u, UBoundaries_u);
		alglib::minlpsetlc2dense(state, D1, LBoundaries_D, UBoundaries_D);
		alglib::minlpsetscale(state, Scale);

		alglib::minlpoptimize(state);
		alglib::minlpresults(state, x, rep);
		
		if (rep.terminationtype <= 4 && rep.terminationtype >= 1) cout << "\nLP OK\n";
		else cout << "\nLP NOT OK\n";
	}


}
