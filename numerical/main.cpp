#include <iostream>
#include <vector>
#include "matrix_class.h"
using namespace std;
int main() {
	Matrix a = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
	Matrix b = {2, 3, 9};
	cout << b << '\n\n';
	//ÐÀÇÎÁÐÀÒÜÑß Ñ ÊÎÍÑÒÐÓÊÒÎÐÎÌ ÊÎÏÈÐÎÂÀÍÈß
	b = b.transpose();
	cout << b(0, 0);
	//cout << b;


	//Matrix ans = GaussSlau(a, b);
	//cout << ans;
	

}