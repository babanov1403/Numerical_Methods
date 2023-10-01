#pragma once
#include <vector>
#include "matrix_class.h"
#include "functions.h"
#include <iostream>


namespace biv {
	class Integrale {
	public:
		int a;
		int b;
		double (*f)(double x) = func;

	};
}
//Interpolation Square Formula?

