#pragma once
#ifndef SPLINE_H
#define SPLINE_H
#include <eigen3/Eigen/Core>

namespace spline
{
class Quintic
{
public:
	// define data type of quintic spline coefficients, Nx6 matrix
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Coef;
	// enumerate three types of boundary conditions
	enum class BC
	{
		Even,
		Odd,
		Mix
	};

public:
	Quintic &operator=(const Quintic &sp);
	Quintic(){};
	Quintic(const Quintic &sp) { *this = sp; }
	// ========== getter ==========
	const Coef &x() const { return _x; };					// x-coordinate (or r)
	const Coef &y() const { return _y; };					// y-coordinate (or z)
	const Eigen::VectorXd &h() const { return _h; };		// chord length
	const Eigen::MatrixX2d &node() const { return _node; }; //spline knots (2-vector)
	// ========== setter ==========
	void node(const Eigen::MatrixX2d &xy); // set node from input matrix
	// compute x- and y-components of quintic spline
	void x(BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	void y(BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	// ========== API ==========
	// compute arc length of i-th local quintic spline
	double localArc(int i, double t = 1.0, int nqd = 20) const;
	// find internal coordiante t of i-th local spline corresponding given an arclength value
	double arc2t(int i, double arc, double eps = 5.e-14, int nqd = 20) const;
	// i-th order differentiation of spline
	const Eigen::Vector3d d(const Coef &x, int i, double t) const;
	// output spline coefficient to file
	void printSpline(const std::string &name);

private:
	// ========== data member ==========
	Coef _x, _y;			// x and y spline coef
	Eigen::VectorXd _h;		// chord
	Eigen::MatrixX2d _node; // node_ (input)

	void h(const Eigen::MatrixX2d &node);
	void setComponent(int i, Coef &x);
	void computeCoef(Coef &x, BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	static const double qd_GL_x20[20]; // Gauss-Legendre quadrature abscissa
	static const double qd_GL_w20[20]; // Gauss-Legendre quadrature weights
};
} // namespace spline

#endif