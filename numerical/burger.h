#pragma once
#include "polynom_class.h"
#include "matrix_class.h"
#include <iostream>
#include "burger_prot.h"
#include "integrate_methods.h"
using namespace std;

namespace biv{

vector<double> makeNodes_burger(int n, double a, double b) {
	vector<double> nodes(n+1);
	double h = (b - a) / n;
	nodes[0] = a;
	for (int i = 1; i <= n; i++)
		nodes[i] = nodes[i - 1] + h;
	return nodes;
}


//works 100% well
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
const int nodes_count = 25;
Matrix<double> A = { {-1, 0}, {0, 3} };
Polynom<Matrix<double>> expRowA = buildRowexp_A(A);
Polynom<Matrix<double>> expA_inv = buildRowexp_A(A * (-1));
Matrix<double> B = Matrix<double>{ 1, 1 }.transpose();
Matrix<double> x0 = Matrix<double>{ 2, 2 }.transpose();
double t0 = 0, t1 = 2;
Matrix<double> expA_t0 = expRowA.value(t0);
Matrix<double> expA_t1 = expRowA.value(t1);
Matrix<double> expA_t1_inv = expA_inv.value(t1);
Matrix<double> expA_t0_inv = expA_inv.value(t0);
size_t N = 2;
Matrix<double> C = Matrix<double>{ 2, 3 }.transpose();
Matrix<double> H = {1, 3};
Matrix<double> g = Matrix<double>{ -1}.transpose();

}

namespace biv {

//две функции - перва€ считает че стоит под интегралом у C, втора€ у D(в точке)

double computeMoment_burger(int k, double a, double b) {
	if (k == -1) {
		cout << "k == -1\n";
		throw;
	}
	double ans_scalar = my_pow(b, k + 1)*1.0 / (k + 1) - my_pow(a, k + 1)*1.0 / (k + 1);
	return ans_scalar;
}
/*
*  огда € считаю Ch => т.е у мен€ из произв матриц получаетс€ скал€р, дл€ которого у мен€ будут
* уже готовые Ai - тоже скал€ры
* потому что в первом случае весова€ функци€ p(x) == 1, а во втором(Dh) она равна p(x) == 1
* 
*/

template <typename Number>
Matrix<Number> makeISF_CD(const vector<Number>& nodes) {
	Matrix<Number> Mu(nodes.size(), 1, 0);
	double a = nodes[0], b = nodes.back();
	for (int i = 0; i < nodes.size(); i++)
		Mu(i, 0) = computeMoment_burger(i, a, b);
	Matrix<Number> A(nodes.size(), nodes.size());
	for (int i = 0; i < nodes.size(); i++)
		for (int j = 0; j < nodes.size(); j++)
			A(i, j) = my_pow(nodes[j], i); 
	Matrix<Number> Outp = GaussSlau(A, Mu);
	return Outp;
}



template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& x0, double t) {
	Matrix<Number> x = expRowA.value(t);
	x = x*expA_t0_inv; 
	x = x * x0;
	return x;
}
template <typename Number>
Matrix<Number> compute_g0(const Matrix<Number>& g, const Matrix<Number>& H, const Matrix<Number>& x0) {
	Matrix<Number> g0 = g;
	Matrix<Number> helper = H;
	helper = helper*compute_x0(x0, t1);
	g0 = g0 + helper * (-1);
	return g0;
}
template <typename Number>
Matrix<Number> compute_G(const Matrix<Number>& H, double t) {
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
		C(0, k) = computeCh(i, i + h, B);
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
	, const Matrix<Number>& B
	, const Matrix<Number>& H) {
	Matrix<Number> res = compute_G(H, a) * B;
	return res;
}

template <typename Number> 
Number func_Ch(
	double a
	, const Matrix<Number>& B) {
	Matrix<Number> res = C.transpose();
	res = res * expA_t1;
	res = res * expA_inv.value(a);
	res = res * B;
	return res[0][0];
}
template <typename Number>
double computeCh(
	double a
	, double b
	, const Matrix<Number>& B) {
	vector<double> nd = makeNodes_burger(nodes_count, a, b);
	Matrix<Number> Aj = makeISF_CD(nd);
	int n = Aj.n;
	double integrale= 0;
	vector<double> sl(n);
	for (int i = 0; i < n; i++)
		sl[i] = Aj(i, 0) * func_Ch(nd[i], B);
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
	vector<double> nd = makeNodes_burger(nodes_count, a, b);
	Matrix<Number> Aj = makeISF_CD(nd);
	int n = Aj.n;
	Matrix<Number> integrale = func_Dh(nd[0], B, H) * Aj[0][0];
	for (int i = 1; i < n; i++) 
		integrale += func_Dh(nd[i], B, H) * Aj[i][0];// Mx1
	return integrale;
}
//x(Tk+1) = Y(Tk+1)*Y^-1(Xk)*x0 + Y(Tk+1) * integrale(xk -> xk+1, Y^-1(t)*b dt)

template <typename Number>
Matrix<Number> func_Cauchy(
	double a
	, const Matrix<Number>& B
	, double u) {
	Matrix<Number> res = expA_inv.value(a);
	res = res * B;
	res = res * u;
	return res;
}

template <typename Number>
Matrix<Number> ComputeIntegraleCauchy(double a
	, double b
	, const Matrix<Number>& B
	, double u) {
	vector<double> nd = makeNodes_burger(nodes_count, a, b);
	Matrix<Number> Aj = makeISF_CD(nd);
	int n = Aj.n;
	Matrix<Number> integrale = func_Cauchy(nd[0], B, u) * Aj[0][0];
	for (int i = 1; i < n; i++)
		integrale += func_Cauchy(nd[i], B, u) * Aj[i][0];
	return integrale;
}
template <typename Number>
Matrix<Number> computePathWGivenU(const vector<Number>& u) {
	double h = (t1 - t0)*1.0 / N;
	Matrix<Number> curr = x0;
	int k = 0;
	for (double i = h; i <= t1; i += h, k++) {
		curr = expRowA.value(i) * (expA_inv.value(i - h) * curr + ComputeIntegraleCauchy(i - h, i, B, u[k]));
		cout << "Current:\n";
		cout << curr;
	}
		
	return curr;
}
/*
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
}*/
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
