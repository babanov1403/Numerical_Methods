#pragma once
#include "matrix_class.h"
using namespace biv;
namespace biv{
template <typename Number>
Matrix<Number> compute_G(const Matrix<Number>& A, const Matrix<Number>& H, double t);
template <typename Number>
Matrix<Number> compute_x0(const Matrix<Number>& A, Matrix<Number>& x0, double t, double t0);
template <typename Number>
Matrix<Number> compute_g0(const Matrix<Number>& g, const Matrix<Number>& H, const Matrix<Number>& x0);
template <typename Number>
Matrix<Number> buildC(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0);
template <typename Number>
Matrix<Number> buildD(
	double t0
	, double t1
	, size_t N
	, const Matrix<Number>& A
	, const Matrix<Number>& b
	, const Matrix<Number>& H
	, const Matrix<Number>& g
	, const Matrix<Number>& x0);
template <typename Number>
double computeCh(
	double a
	, double b
	, Matrix<Number> A
	, const Matrix<Number>& B);
template <typename Number>
Matrix<Number> computeDh(
	double a
	, double b
	, double t1
	, const Matrix<Number>& B
	, const Matrix<Number>& A
	, const Matrix<Number>& H);
}
