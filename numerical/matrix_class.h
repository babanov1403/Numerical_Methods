#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;
namespace biv {
	class Matrix;
	ostream& operator<<(ostream& os, const Matrix& a);

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
			: v{ init }
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
		num& operator()(int i, int j) {
			if (i >= this->n || j >= this->m) {
				cout << "Seg Fault\n";
				throw;
			}
			return this->v[i][j];
		}

		const num& operator()(int i, int j) const {
			if (i >= this->n || j >= this->m) {
				cout << "Seg Fault\n";
				throw;
			}
			return this->v[i][j];
		}
		//Multiply by number
		Matrix operator*(double x) {
			Matrix m = this->v;
			for (auto& i : m.v)
				for (auto& j : i) j *= x;
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
		Matrix transpose() {
			Matrix& a = *this;
			Matrix b(a.m, a.n);
			for (int i = 0; i < a.m; i++)
				for (int j = 0; j < a.n; j++)
					b.v[i][j] = a.v[j][i];
			return b;
		}
	};
	ostream& operator<<(ostream& os, const Matrix& a) {
		for (auto i : a.v) {
			for (double j : i)
				os << j << " ";
			os << '\n';
		}
		return os;
	}

	Matrix GaussSlau(const Matrix& A, const Matrix& b) {
		size_t n = b.n;
		if (A.m != A.n) {
			cout << "Matrix A is not a square!\n";
			throw;
		}
		Matrix W(n, n + 1);
		for (int i = 0; i < n; i++) {
			copy((A.v[i]).begin(), (A.v[i]).end(), (W.v[i]).begin());
			W.v[i][n] = b.v[i][0];
		}
		for (size_t j = 0; j < n; j++) {
			int maxx = j;
			for (size_t i = j + 1; i < n; i++)
				if (abs(W(maxx, j)) > abs(W(i, j))) {
					swap(W.v[maxx], W.v[i]);
					maxx = i;
				}
			double denom = W(j, j);
			for (size_t i = 0; i < n + 1; i++) {
				W(j, i) = W(j, i) * 1.0 / denom;
			}
			for (size_t i = j + 1; i < n; i++) {
				denom = W(i, j);
				for (size_t q = 0; q < n + 1; q++) {
					W(i, q) = W(i, q) - (W(j, q) * denom);
				}
			}
		}
		// inverse
		for (int j = n - 1; j > 0; j--)
		{
			//for (int dd = 0; dd < n; dd++)
			//	W.v[j][dd] /= W.v[j][j];

			for (int i = j - 1; i >= 0; i--) {
				double denom = W(i, j);
				for (int dd = 0; dd < n + 1; dd++) //CAN GET MORE EFFIENCY
					W(i, dd) = W(i, dd) - (W(j, dd) * denom);
			}
		}
		Matrix res(n, 1);
		for (size_t i = 0; i < n; i++)
			res.v[i][0] = W.v[i][n];

		return res;
	}
}