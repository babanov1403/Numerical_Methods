#pragma once
#include "polynom_class.h"
#include "matrix_class.h"
#include <iostream>
#include "integrate_methods.h"
#include "burger_prot.h"
using namespace std;

namespace biv{
template <typename Number>
Polynom<Matrix<Number>> buildRowexp_A(const Matrix<Number>& A, int k = 25) {
	if (A.n != A.m) {
		std::cerr << "what do you think i will do with this matrix?\n";
		throw;
	}
	Matrix<Number> XD(A.n, A.n, 0);
	for (int i = 0; i < A.n; i++) {
		XD(i, i) = 1;
	}
	Polynom<Matrix<Number>> res(k);
	double curr = 1;
	res[0] = XD*curr;
	for (int i = 1; i < k; i++) {
		curr *= 1.0 / i;
		XD = XD*A;
		res[i] = XD * curr;
	}
	return res;
}
Matrix<double> A = { {1, 0}, {0, -2} };
static Polynom<Matrix<double>> expRowA = buildRowexp_A(A);
Matrix<double> B = Matrix<double>{ 1, 2 }.transpose();
Matrix<double> x0 = Matrix<double>{ 1, 2 }.transpose();
double t0 = 0, t1 = 10;
static Matrix<double> expA_t0 = expRowA.value(t0);
static Matrix<double> expA_t1 = expRowA.value(t1);
static Matrix<double> expA_t1_inv = expA_t1.inverse();
size_t N = 2;
Matrix<double> C = Matrix<double>{ 4, 1 }.transpose();
Matrix<double> H = { {2, 1}, {3, 4} };
Matrix<double> g = Matrix<double>{ -1, 1 }.transpose();
}

namespace biv {

//две функции - первая считает че стоит под интегралом у C, вторая у D(в точке)



template <typename Number>
Matrix<Number> computeMoment(int k, int matrix_size, double a, double b) {
	if (k == -1) {
		cout << "k == -1\n";
		throw;
	}
	Matrix<Number> ans(matrix_size);
	double ans_scalar = my_pow(b, k + 1)*1.0 / (k + 1) - my_pow(a, k + 1)*1.0 / (k + 1);
	for (int i = 0; i < matrix_size; i++)
		ans[i][i] = ans_scalar;
	return ans;
}
/*
* Когда я считаю Ch => т.е у меня из произв матриц получается скаляр, для которого у меня будут
* уже готовые Ai - тоже скаляры
* потому что в первом случае весовая функция p(x) == 1, а во втором(Dh) она равна p(x) == 1
*/

template <typename Number>
Matrix<Number> makeISF_CD(const vector<Number>& nodes) {
	Matrix<Number> Mu(nodes.size(), 1, 0);
	double a = nodes[0], b = nodes.back();
	for (int i = 0; i < nodes.size(); i++)
		Mu(i, 0) = computeMoment(i, 1, a, b);
	Matrix<Number> A(nodes.size(), nodes.size());
	for (int i = 0; i < nodes.size(); i++)
		for (int j = 0; j < nodes.size(); j++)
			A(i, j) = my_pow(nodes[j], i); 
	Matrix<Number> Outp = GaussSlau(A, Mu);
	return Outp;
}



template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& A, Matrix<Number>& x0, double t, double t0) {
	Matrix<Number> x = expRowA.value(t);
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
	Matrix<Number> G = H;
	G = G * expA_t1;
	G = G * (expRowA.value(t)).inverse();
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
Matrix<Number> func_Dh(
	double a
	, double t1
	, const Matrix<Number>& B
	, const Matrix<Number>& A
	, const Matrix<Number>& H) {
	Matrix<Number> res = compute_G(A, H, a) * B;
	return res;
}

template <typename Number> 
Number func_Ch(
	double a
	, Matrix<Number> A
	, const Matrix<Number>& B) {
	Matrix<Number> res = C.transpose();
	res = res * expRowA.value(a);
	res = res * expA_t1_inv;
	res = res * B;
	return res;
}

//доделать считатели и продебажить
template <typename Number>
double computeCh(double a
	, double b
	, Matrix<Number> A
	, const Matrix<Number>& B) {
	int n = Aj.n;
	double integrale;
	vector<double> sl(n);
	for (int i = 0; i < n; i++)
		sl[i] = Aj(i, 0) * func_Ch(nodes[i]);
	sort(sl.begin(), sl.end());
	for (auto i : sl) integrale += i;
	return integrale;
}
template <typename Number>
Matrix<double> computeDh(double a
	, double b
	, double t1
	, const Matrix<Number>&B
	, const Matrix<Number>&A
	, const Matrix<Number>&H) {
	int n = Aj.n;
	Matrix<double> integrale(;
	vector<double> sl(n);
	for (int i = 0; i < n; i++) {
		sl[i] = Aj(i, 0) * func_Dh(nodes[i]);// xj = tj + a(a a + b / 2 b)
	}
	sort(sl.begin(), sl.end());
	for (auto i : sl) integrale += i;
	return integrale;
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
	res = res * expRowA.value(a);
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
		curr = curr * expRowA.value(ksi).inverse();
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
