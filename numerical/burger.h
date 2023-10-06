#pragma once
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"


template <typename Number>
Matrix<Number> buildRow(const Matrix<Number>& A, double t, int k = 100) {
	Matrix<Number> T(1, 1, t);
	Matrix<Number> Res(A.n, A.n, 0);
	for (int q = 0; q <= k; q++) {
		Res += A ^ q * T ^ q * (1.0 / fact(q));
	}
}