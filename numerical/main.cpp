#include <iostream>
#include <vector>
#include "integrate_methods.h"
#include "matrix_class.h"
#include <cmath>
using namespace std;
using namespace biv;
int main() {
	Matrix A = { {1, 1, 1}, {1, 2, 3}, {1,4,9} };
	Matrix b = { 2 * sqrt(2), 4 / 3 * sqrt(2), 8 / 5 * sqrt(2) };
	b = b.transpose();
	Matrix ans = GaussSlau(A, b);
	cout << ans;
	

}