#pragma once
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"

namespace biv {

template <typename Number>
Matrix<Number> buildRowexpA(const Matrix<Number>&, double, int);
template <typename Number>
Matrix<Number> compute_G(const Matrix<Number>& A, const Matrix<Number>& H, double t, double t1);
template <typename Number>
Matrix<Number> computeMatrixMulti(const vector<Matrix<Number>>& v, double ksi);
template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& A, Matrix<Number>& x0, double t, double t0);
template <typename Number>
Matrix<Number> compute_g0(const Matrix<Number>& g, const Matrix<Number>& H, const Matrix<Number>& x0);
template <typename Number>
Matrix<Number> buildC(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0);
template <typename Number>
Matrix<Number> buildD(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0);
template <typename Number>
double computeCh(
	double a
	, double b
	, Matrix<Number> A
	, const Matrix<Number>& B);
template <typename Number>
Matrix<Number> computeDh(
	double a
	, double b
	, double t1
	, const Matrix<Number>& B
	, const Matrix<Number>& A
	, const Matrix<Number>& H);


template <typename Number>
Matrix<Number> buildRowexpA(const Matrix<Number>& A, double t, int k = 30) {
	Matrix<Number> T(1, 1, t);
	Matrix<Number> Res(A.n, A.n, 0);
	Matrix<Number> CurrA(A.n, A.n, 0);
	double currnumber = 1;
	//https://habr.com/ru/articles/239303/

	for (int i = 0; i < A.n; i++) CurrA(i, i) = 1;
	Res += CurrA;
	for (double q = 1; q <= k; q++){//might be better
		CurrA = CurrA*A;
		currnumber *= 1.0 / q;
		Res += CurrA*currnumber;
	}
	Res = Res^t;
	return Res;
}
template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& A, Matrix<Number>& x0, double t, double t0) {
	Matrix<Number> x = buildRowexpA(A, t);
	x = x*buildRowexpA(A, t0).inverse(); 
	x = x * x0;
	return x;
}
template <typename Number>
Matrix<Number> compute_g0(const Matrix<Number>& g, const Matrix<Number>& H, const Matrix<Number>& x0) {
	Matrix<Number> g0 = g;
	Matrix<Number> helper = H;
	helper = helper*x0;
	g0 = g0 + helper * (-1);
	return g0;
}
template <typename Number>
Matrix<Number> compute_G(const Matrix<Number>& A, const Matrix<Number>& H, double t, double t1) {
	Matrix<Number> G = H;
	G = G * buildRowexpA(A, t1);
	G = G * (buildRowexpA(A, t).inverse());
	return G;
}
template <typename Number>
Matrix<Number> buildC(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0) {
	if (t1 <= t0) {
		cout << "bad build C\n";
		throw;
	}
	Matrix<Number> C(1, N);
	double h = (t1 - t0) / N;
	int k = 0;
	for (double i = t0; i <= t1-h; i += h, k++) {
		C(0, k) = computeCh(i, i + h, A, b);
	}
	return C;
}
template <typename Number>
Matrix<Number> buildD(
	  double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0) {
	if (t1 <= t0) {
		cout <<"bad build D\n";
		throw;
	}
	double h = (t1 - t0) / N;
	Matrix<Number> D(g.n, N, 0);
	int q = 0;
	for (double i = t0; i <= t1-h; i += h, q++) {
		Matrix<Number> curr = computeDh(i, i + h, t1, b, A, H);
		for (int k = 0; k < g.n; k++) {
			D(k, q) = curr(k, 0);
		}
	}
	return D;
}
template <typename Number>
Matrix<Number> computeDh(
	  double a
	, double b
	, double t1
	, const Matrix<Number>& B
	, const Matrix<Number>& A
	, const Matrix<Number>& H) {
	if (a >= b) {
		cout << "bad Dh\n";
		throw;
	}
	double h = (b - a) / 1e3;
	Matrix<Number> res = compute_G(A, H, a, t1);
	res = res * B;
	res *= h;

	Matrix<Number> curr = res;
	cout << "ksi:\n";
	for (double ksi = a + h; ksi <= b; ksi += h) {
		cout << ksi << '\n';
		curr = compute_G(A, H, ksi, t1) * B; //а почему дискретно то емае
		cout << curr;
		curr *= h;
		res += curr;
	}
	cout << a << " " << b << '\n';
	cout << res;
	return res;
}
template <typename Number>
double computeCh(
	  double a
	, double b
	, Matrix<Number> A
	, const Matrix<Number>& B) {
	if (a >= b) {
		cout << "bad Ch\n";
		throw;
	}
	double h = (b - a) / 1e3;
	A = A * (-1);
	A = A.transpose();
	Matrix<Number> res = buildRowexpA(A, a);
	res = res * B;
	res *= h;
	Matrix<Number> curr = res;
	for (double ksi = a+h; ksi <= b; ksi += h) {
		curr = buildRowexpA(A, ksi)*B;
		curr *= h;
		res += curr;
	}
	double ans = res(0, 0);
	return ans;
}
/*template <typename Number>
Matrix<Number> RhimannIntegrale(double x0, double x1, int N = 1e5, Matrix<Number> (*MatrixFunc)(double ksi)) {
	double h = (x1 - x0) / N;
	Matrix<Number> res = MatrixFunc(x0);
	Matrix<Number> curr = res;
	for (double ksi = x0+h; ksi <= x1; ksi += h) {
		curr = MatrixFunc(ksi);
		curr *= h;
		res += h;
	}
	return res;
}*/
}
