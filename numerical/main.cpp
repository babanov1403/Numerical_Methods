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
	
	


}