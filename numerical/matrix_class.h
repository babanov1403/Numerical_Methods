#pragma once
#include <iostream>
#include <vector>
using namespace std;
class Matrix {
	using num = double;
public:
	vector<vector<num>> v;
	size_t n;
	size_t m;

	Matrix() : n(0), m(0) {}
	Matrix(const vector<vector<num>>& v)
		: v(v) 
		, n(v.size())
		, m(v[0].size()) {}
	Matrix(const vector<num>& v) : n(1), m(v.size()) {
		vector<vector<num>> copy(1, vector<num>(v.size()));
		copy[0] = v;
		this->v = copy;
	}
	Matrix(const initializer_list<initializer_list<num>>& init) {
		n = init.size();
		m = (*init.begin()).size();
		for (auto data : init) {
			v.push_back(data);
		}
	}
	Matrix(const initializer_list<num>& init) 
		: v{init}
		, n(1)
		, m(init.size()) {}
	Matrix(size_t n, size_t m)
		: n(n)
		, m(m) {
		vector<vector<num>> q(n, vector<num>(m));
		v = q;
	}
	//copy constr
	Matrix(const Matrix& a) {
		n = a.n;
		m = a.m;
		v = a.v;
	}

	int operator()(int i, int j) {
		return this->v[i][j];
	}
	//Multiply by number
	Matrix operator*(double x) {
		Matrix m = this->v;
		for (auto &i : m.v)
			for (auto &j : i) j *= x;
		return m;
	}
	Matrix& operator*=(double x) {
		Matrix m = this->v;
		for (auto& i : this->v)
			for (auto& j : i) j *= x;
		return m;
	}
	Matrix& operator=(Matrix a) {
		n = a.n;
		m = a.m;
		v = a.v;
		return a;
	}

	//Matrix Multipl
	Matrix operator*(Matrix b) {
		if (this->m != b.n) {
			cout << "Bad dimensions in multiplying\n";
			return b;
		}
		double curr_sum = 0;
		Matrix a = this->v;
		Matrix res(a.n, b.m);
		for (int i = 0; i < a.n; i++) {
			for (int j = 0; j < a.n; j++) {//i, j position in new matrix
				res.v[i][j] = 0;
				for (int q = 0; q < a.m; q++)
					res.v[i][j] += a.v[i][q] * b.v[q][j];
			}
		}
		return res;
	}
	/*
	num& operator[](const int i, const int j) {
		return this->v[i];
	}*/
	Matrix transpose() {
		Matrix& a = *this;
		cout << a;
		Matrix b(a.m, a.n);
		//cout << a.m << " " << a.n << '\n';
		for (int i = 0; i < a.m; i++)
			for (int j = 0; j < a.n; j++)
				b.v[i][j] = a.v[j][i];
		return b;
	}
};

ostream& operator<<(ostream& os, const Matrix& a) {
	for (auto i : a.v) {
		for (auto j : i)
			os << j << " ";
		os << '\n';
	}
	return os;
}
/*
Matrix GaussSlau(const Matrix& A, const Matrix& b) {
	size_t n = b.m;
	cout << n << '\n';
	Matrix W(n, n + 1);
	for (size_t i = 0; i < n; i++) {
		cout << i << '\n';
		copy((A.v[i]).begin(), (A.v[i]).end(), W.v[i]);
		W.v[i][n] = b.v[i][0];
	}
	for (size_t j = 0; j < n; j++) {
		size_t maxx = j;
		for (size_t i = j + 1; i < n; i++)
			maxx = abs(W.v[i][j]) > abs(W.v[i][maxx]) ? j : maxx;

		std::swap(W.v[maxx], W.v[j]);
		for (int dd = 0; dd < W.m; dd++)
			W.v[j][dd] /= W.v[j][j];

		for (size_t i = j + 1; i < n; i++)
			for (size_t dd = 0; dd < n; dd++)
				W.v[i][dd] -= (W.v[j][dd] * W.v[i][j]);
	}
	// inverse
	for (size_t j = n - 1; j > 0; j--)
	{
		for (int dd = 0; dd < n; dd++)
			W.v[j][dd] /= W.v[j][j];

		for (size_t i = j - 1; i < n; i--)
			for (int dd = 0; dd < n; dd++)
				W.v[i][dd] = W.v[i][dd] - (W.v[j][dd] * W.v[i][j]);
	}
	Matrix res(n, 1);
	for (size_t i = 0; i < n; i++)
		res.v[i][0] = W.v[i][n];

	return res;
}*/