#include <iostream>
#include <iomanip>
#include <vector>
#include "integrate_methods.h"
#include "matrix_class.h"
#include <cmath>
#include "burger.h"

using namespace biv;
using namespace std;
vector<double> makeNodes(int n, short int type = 0) {
	if (type == 0) {//RAVNOMERNO
		vector<double> nodes(n);
		double h = (var_b::b - var_b::a) / n;
		nodes[0] = var_b::a;
		for (int i = 1; i < n; i++) 
			nodes[i] = nodes[i - 1] + h;
		return nodes;
	}
}
int main() {
	/*vector<double> nodes = makeNodes(27);
	double integrale = computeIntegral(makeISF(nodes), nodes);
	cout << setprecision(11) << integrale << " ";*/
	
	//cout << "g0:\n";
	//Matrix<double> g0 = compute_g0(g, H, compute_x0(A, x0, t1, t0));
	//cout << g0;
	//cout << "G:\n";
	//Matrix<double> G = compute_G(A, H, t0, t1);
	//cout << G;
	Matrix<double> C_res = buildC(t0, t1, N, A, B, H, g, x0);
	cout << "C:\n";
	cout << C_res;
	//Matrix<double> D_res = buildD(t0, t1, N, A, B, H, g, x0);
	//cout << "D:\n";
	//cout << D_res;
	

	

	
	
}