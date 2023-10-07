#pragma once
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"

namespace biv {

	template <typename Number>
	Matrix<Number> buildRowexpA(const Matrix<Number>&, double, int);
	template <typename Number>
	Matrix<Number> buildG();
	template <typename Number>
	Matrix<Number> computeMatrixMulti(const vector<Matrix<Number>>& v, double ksi);

	template <typename Number>
	Matrix<Number> buildRowexpA(const Matrix<Number>& A, double t, int k = 200) {
		Matrix<Number> T(1, 1, t);
		Matrix<Number> Res(A.n, A.n, 0);
		Matrix<Number> CurrA(A.n, A.n, 0);
		double currnumber = 1;
		//https://habr.com/ru/articles/239303/

		for (int i = 0; i < A.n; i++) CurrA(i, i) = 1;
		for (double q = 1; q <= k; q++) {//might be better
			cout << Res;
			CurrA *= A;
			currnumber *= 1.0 / q;
			cout << CurrA;
			cout << currnumber << "\n";
			Res += CurrA*currnumber;
		}
		cout << Res;
		Res = Res^(t);
		return Res;
	}
	template <typename Number>
	Matrix<Number> compute_x0(const Matrix<Number> A, const Matrix<Number> x0, double t, double t0) {
		Matrix<Number> x = buildRowexpA(A, t);
		cout << x;
		cout << "x^^^^^^^^^\n";
		cout << pow(2.718281828, -20) << "\n";
		//x *= buildRowexpA(A, t0).inverse();
		//x *= x0;
		return x;
	}
	template <typename Number>
	Matrix<Number> compute_g0(const Matrix<Number>& g, const Matrix<Number>& H, const Matrix<Number>& x0) {
		Matrix<Number> g0 = g;
		Matrix<Number> helper = H;
		helper *= x0;
		g0 = g0 + helper * (-1);
		return g0;
	}
}
