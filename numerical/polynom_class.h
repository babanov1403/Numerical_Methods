#pragma once
#include <cmath>
#include <iostream>

template <typename Number>
class Polynom {
public:
	int n;
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
	//
	Number& operator[](int i) {
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
	Polynom& operator*(const Polynom& p) {
		if (p.n > 1) {
			cout << "wait for FFT\n";
			throw;
		}
		else if (p.n == 1) {
			//ÇÀÊÎÍ×ÈË ÇÄÅÑÜ
		}
	}


};