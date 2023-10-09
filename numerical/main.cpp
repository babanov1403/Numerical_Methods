#include <iostream>
#include <iomanip>
#include <vector>
#include "integrate_methods.h"
#include "matrix_class.h"
#include <cmath>
#include "burger.h"

using namespace std;
using namespace biv;
vector<double> makeNodes(int n, short int type = 0) {
	if (type == 0) {//RAVNOMERNO
		vector<double> nodes(n);
		double h = (b - a) / n;
		nodes[0] = a;
		for (int i = 1; i < n; i++) 
			nodes[i] = nodes[i - 1] + h;
		return nodes;
	}
}
int main() {
	/*vector<double> nodes = makeNodes(27);
	double integrale = computeIntegral(makeISF(nodes), nodes);
	cout << setprecision(11) << integrale << " ";*/
	/*Matrix<double> GGG = makeSFGauss<double>(3);
	cout << GGG;*/
	Matrix<double> A = { {1, 0}, {0, -2} };
	Matrix<double> b = Matrix<double>{ 1, 2 }.transpose();
	Matrix<double> x0 = Matrix<double>{ 1, 2 }.transpose();
	cout << "x0:\n";
	cout << x0;
	double t0 = 0, t1 = 10;
	double h = 1;
	Matrix<double> C = Matrix<double>{ 4, 1 }.transpose();
	Matrix<double> H = { {2, 1}, {3, 4} };
	Matrix<double> g = Matrix<double>{ -1, 1 }.transpose();
	Matrix<double> g0 = compute_g0(g, H, compute_x0(A, x0, t1, t0));
	cout << "g0:\n";
	cout << g0;

	/*Matrix<double> Asopr = A * (-1);
	Asopr = Asopr.transpose();*/



	//Õ¿ƒŒ ¡”ƒ≈“ œŒ—“–Œ»“‹ Ã¿“–»÷” G(t) » Õ¿…“» ‘”Õƒ¿Ã≈Õ“¿À‹Õ”ﬁ ƒÀﬂ -A^T » œŒ—À≈ ›“Œ√Œ œŒ—◊»“¿“‹ »Õ“≈√–¿À€
	//» ¬—≈ √√ ¬€¬≈—“»  –¿—»¬Œ ¬≈ “Œ– U(u1, u2, u3, u4....., un) Ë Â˘Â ÔÓÒ˜ËÚ‡Ú¸ g0 => Ò‰ÂÎ‡ÌÓ
	
}