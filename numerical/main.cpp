#include <iostream>
#include <vector>
#include "matrix_class.h"
using namespace std;
int main() {
	Matrix a = { {1, -2, 3}, {3, 1, -1}, {2, 5, 2} };
	Matrix b = {2, 3, 9};
	cout << a;
	cout << b;
	b = b.transpose();

	Matrix ans = GaussSlau(a, b);

	cout << ans;
	

}