#pragma once
#include "matrix_class.h"
#include <cmath>
#include <string>
#include <iostream>

namespace biv {

double my_pow(double x, int k) {
	double res = 1;
	while (k > 0) {
		if ((k & 1)) res = res * x;
		x *= x;
		k >>= 1;
	}
	return res;
}

template <typename Number>
class Polynom {
public:
	size_t n;
	vector<Number> v;

	Polynom(int n) : n(n) {
		if (n < 0) {
			cout << "Bad degree\n";
			throw;
		}
		vector<Number> a(n);
		v = a;
	}
	Polynom(const vector<Number>& v) : v(v) {}
	Polynom(const std::initializer_list<Number>& ll) {
		v = ll;
		n = ll.size();
	}
	Polynom(const Polynom& p) {
		v = p.v;
		n = p.n;
	}

	Number& operator[](int i) {
		if (i >= n || i < 0) {
			cout << "Seg fault poly\n";
			throw;
		}
		return &v[i];
	}
	const Number& operator[](int i) const {
		if (i >= n || i < 0) {
			cout << "Seg fault poly\n";
			throw;
		}
		return v[i];
	}
	Polynom& operator*(Number x) {
		for (auto& q : v) q *= x;
		return *this;
	}
	/// <summary>
	/// ??????????????
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	Polynom& operator*(const Polynom& p) {
		if (p.n > 1) {
			cout << "wait for FFT\n";
			throw;
		}
		else if (p.n == 1) {
			//�������� �����
		}
		return {};
	}

	Number value(double x) {//O(nlogn) => i can do this by O(n)

		Number ans = 0;
		if (abs(x) <= 1) {
			for (int i = this->n - 1; i >= 0; i--)
				ans = ans + this[i] * my_pow(x, i);
			return ans;
		}
		else {
			for (int i = 0; i < this->n; i++)
				ans = ans + this[i] * my_pow(x, i);
			return ans;
		}
	}
};


template <typename Number>
class Polynom<Matrix<Number>> {
public:
	size_t n;
	vector<Matrix<Number>> v;

	Polynom(int n) : n(n) {
		if (n < 0) {
			cout << "Bad degree\n";
			throw;
		}
		vector<Matrix<Number>> a(n);
		v = a;
	}
	Polynom(const vector<Matrix<Number>>& v) : v(v) {}
	Polynom(const std::initializer_list<Matrix<Number>>& ll) {
		v = ll;
		n = ll.size();
	}
	Polynom(const Polynom& p) {
		v = p.v;
		n = p.n;
	}

	Matrix<Number>& operator[](int i) {
		if (i >= n || i < 0) {
			cout << "Seg fault poly\n";
			throw;
		}
		return v[i];
	}
	const Matrix<Number>& operator[](int i) const {
		if (i >= n || i < 0) {
			cout << "Seg fault poly\n";
			throw;
		}
		return v[i];
	}
	Polynom& operator*(Number x) {
		for (auto& q : v) q *= x;
		return *this;
	}
	/// <summary>
	/// ??????????????
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	Polynom& operator*(const Polynom& p) {
		if (p.n > 1) {
			cout << "wait for FFT\n";
			throw;
		}
		else if (p.n == 1) {
			//??????????
		}
		return {};
	}

	Matrix<Number> value(double x) {
		if (n == 0) {
			cerr << "Empty poly\n";
			throw;
		}
		Matrix<Number> ans(v[0].n, v[0].n, 0);
		Matrix<Number> tmp = ans;
		if (abs(x) <= 1) {
			for (int i = n - 1; i >= 0; i--) {
				tmp = v[i] * my_pow(x, i);
				ans = ans + tmp;
			}
			return ans;
		}
		else {
			int helper = floor(x);
			x /= helper;
			for (int i = 0; i < n; i++) {
				tmp = v[i] * my_pow(x, i);
				ans = ans + tmp;
			}
			ans = ans ^ helper;
			return ans;
		}
	}
};

template <typename Number>
ostream& operator<<(ostream& os, const Polynom<Number>& poly) {
	if (poly.n == 0) {
		os << "Empty poly";
		return os;
	}
	for (int i = poly.n - 1; i >= 0; i--) {
		std::string prefix = "";
		std::string suffix = "";
		if (poly[i] >= 0 && i != poly.n - 1) prefix += '+';
		else if (i != poly.n - 1) prefix += '-';
		os << prefix;
		os << poly[i];
		if (i != 0)
			suffix += "x^" + std::to_string(i);
		os << suffix;
	}
	return os;
}

template <>
ostream& operator<< <Matrix<double>>(ostream& os, const Polynom<Matrix<double>>& poly) {
	if (poly.n == 0) {
		os << "Empty poly";
		return os;
	}
	//??? do smth pretty
	for (int i = poly.n - 1; i >= 0; i--) {
		std::string suffix = "";
		os << poly[i];
		if (i != 0)
			suffix += "x^" + std::to_string(i);
		os << suffix << '\n';
	}
	os << "\n A 4to ti hotel ? ? ? ? \n";
	return os;
}
}