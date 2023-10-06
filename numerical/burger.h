#pragma once
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"

namespace biv {
	static const double t0 = 0;
	static Matrix<double> X0 = Matrix<double>{ 1, 2 }.transpose();

	template <typename Number>
	Matrix<Number> buildRowexpA(const Matrix<Number>& A, double t, int k = 20) {
		std::cout << X0;
		Matrix<Number> T(1, 1, t);
		Matrix<Number> Res(A.n, A.n, 0);
		Matrix<Number> Curr(A.n, A.n, 0);
		//maybe do backwards?
		for (int i = 0; i < A.n; i++) Curr(i, i) = 1;
		for (int q = 1; q <= k; q++) {
			cout << Res;
			cout << '\n';
			Curr *= A;
			Curr *= t / q;
			Res += Curr;
		}
		return Res;
	}
	/*
	*x' = Ax + bu
	*x(t0) = x0
	* [t0, t1]
	* h = (t1 - t0)/N
	* vector<int> U;??? Matrix?
	* C.transpose()*x(t1) -> max
	* H * x(t1) = g
	* 1) Ch(t) = integrale(Z^T*b(t)dt)
	*
	* 
	* 
	*/
}
