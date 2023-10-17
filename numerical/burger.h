#pragma once
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"
#include "burger_prot.h"

namespace biv {
template <typename Number>
Matrix<Number> buildRowexpA(const Matrix<Number>& A, double t, int k = 300) {
	Matrix<Number> T(1, 1, t);
	Matrix<Number> Res(A.n, A.n, 0);
	Matrix<Number> CurrA(A.n, A.n, 0);
	double currnumber = 1;
	//https://habr.com/ru/articles/239303/
	for (int i = 0; i < A.n; i++) CurrA(i, i) = 1;
	Res += CurrA;
	int counter = 0;
	while (t > 10) {//might be better
		t /= 2;
		k++;
	}
	for (double q = 1; q <= k; q++) {//might be better
		CurrA = CurrA * A;
		currnumber *= t / q;
		Res += CurrA * currnumber;
	}
	Res = Res ^ static_cast<int>(pow(2, counter));
	return Res;
}
Matrix<double> A = { {1, 0}, {0, -2} };
Matrix<double> B = Matrix<double>{ 1, 2 }.transpose();
Matrix<double> x0 = Matrix<double>{ 1, 2 }.transpose();
double t0 = 0, t1 = 10;
static Matrix<double> expA_t0 = buildRowexpA(A, t0);
static Matrix<double> expA_t1 = buildRowexpA(A, t1);
static Matrix<double> expA_t1_inv = expA_t1.inverse();
size_t N = 2;
Matrix<double> C = Matrix<double>{ 4, 1 }.transpose();
Matrix<double> H = { {2, 1}, {3, 4} };
Matrix<double> g = Matrix<double>{ -1, 1 }.transpose();
}

namespace biv {

template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& A, Matrix<Number>& x0, double t, double t0) {
	Matrix<Number> x = buildRowexpA(A, t);
	x = x*expA_t0.inverse(); 
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
Matrix<Number> compute_G(const Matrix<Number>& A, const Matrix<Number>& H, double t) {
	//cout << t << "<-G->" << t1 << '\n';
	Matrix<Number> G = H;
	//cout << G;
	//cout << expA_t1;
	G = G * expA_t1;
	//cout << G;
	//cout << "exp(A):\n";
	//cout << buildRowexpA(A, t);
	G = G * (buildRowexpA(A, t).inverse());
	//cout << "G:\n";
	//cout << G;
	return G;
}
template <typename Number>
Matrix<Number> buildC(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& B
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
		C(0, k) = computeCh(i, i + h, A, B);
	}
	return C;
}
template <typename Number>
Matrix<Number> buildD(
	  double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& B
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
		//cout << "FIRST ITER:\n";
		Matrix<Number> curr = computeDh(i, i + h, t1, B, A, H);
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
	//cout << a << "<->" << b << '\n';
	double h = (b - a) / 1e4;
	Matrix<Number> res = compute_G(A, H, a)*B;
	res *= h;
	Matrix<Number> curr = res;
	for (double ksi = a + h; ksi <= b-h; ksi += h) {
		curr = compute_G(A, H, ksi);
		curr = curr * B;
		curr *= h;
		res += curr;
	}
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
	double h = (b - a) / 1e5;
	//cout << A;
	Matrix<Number> res = C.transpose();
	//cout << "DOT:\n";
	//cout << a << '\n';
	res = res * buildRowexpA(A, a);
	res = res * expA_t1_inv;
	res = res * B;
	res *= h;
	Matrix<Number> curr;
	for (double ksi = a+h; ksi < b-h; ksi += h) {
		//cout << ksi << "\n";
		curr = C.transpose();
		//cout << "xC:\n";
		//cout << curr;
		curr = curr * expA_t1;
		//cout << "xe^A(ksi):\n";
		//cout << curr;
		curr = curr * buildRowexpA(A, ksi).inverse();
		//cout << "xe^A-1(t*):\n";
		//cout << curr;
		curr = curr * B;
		//cout << "B:\n";
		//cout << curr;
		curr *= h;
		res += curr;
	}
	cout << res;
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
