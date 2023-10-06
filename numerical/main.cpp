#include <iostream>
#include <vector>
#include "integrate_methods.h"
#include "matrix_class.h"
#include <cmath>

using namespace std;
using namespace biv;
int main() {
	Matrix<double> A = { {1, -1}, {2, 4} };
	A = A.inverse();
}