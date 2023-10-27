#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "integrate_methods.h"
#include "matrix_class.h"
#include "burger.h"
#include "polynom_class.h"



using namespace biv;
using namespace std;
int main() {
	Matrix<double> C_res = buildC(t0, t1, N, A, B, H, g, x0);
	cout << "C:\n";
	cout << C_res;
	Matrix<double> D_res = buildD(t0, t1, N, A, B, H, g, x0);
	cout << "D:\n";
	cout << D_res;
	cout << "g0:\n";
	cout << compute_g0(g, H, x0);


}
//Matrix<double> D_res = buildD(t0, t1, N, A, B, H, g, x0);
	//cout << "D:\n";
	//cout << D_res;
	/*vector<double> nodes = makeNodes(27);
	double integrale = computeIntegral(makeISF(nodes), nodes);
	cout << setprecision(11) << integrale << " ";*/

	//cout << "g0:\n";
	//Matrix<double> g0 = compute_g0(g, H, compute_x0(A, x0, t1, t0));
	//cout << g0;
	//cout << "G:\n";
	//Matrix<double> G = compute_G(A, H, t0, t1);
	//cout << G;