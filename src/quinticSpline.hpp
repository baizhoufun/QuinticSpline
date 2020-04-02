#pragma once
#ifndef SPLINE_H
#define SPLINE_H
#include <eigen3/Eigen/Core>

namespace spline
{
class Quintic
{
private:
	// define data type of quintic spline coefficients, Nx6 matrix
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Coef;
	// ========== data member ==========
	std::vector<Coef> _component;
	Eigen::VectorXd _h;	   // chord
	Eigen::MatrixXd _node; // node_ (input)
	int _dim = 2;

public:
	// enumerate three types of boundary conditions
	enum class BC
	{
		Even,
		Odd,
		Mix
	};
	Quintic &operator=(const Quintic &sp);
	Quintic(){};
	Quintic(const Quintic &sp) { *this = sp; }
	const int &dim() const { return _dim; };
	const std::vector<Coef> &component() const { return _component; };	// x-coordinate (or r)
	const Coef &x() const { return _component[0]; };					// x-coordinate (or r)
	const Coef &y() const { return _component[1]; };					// x-coordinate (or r)
	const Coef &z() const { return _component[2]; };					// x-coordinate (or r)
	const Eigen::VectorXd &h() const { return _h; };					// chord length
	const Eigen::MatrixXd &node() const { return _node; };				//spline knots
	void node(const Eigen::MatrixXd &xy);								// set node from input matrix
	void node(const Eigen::MatrixXd &xy, const Eigen::VectorXd &chord); // set node from input matrix
	void computeComponent(int k, BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	// ========== API ==========
	const Eigen::VectorXd arcCoord() const;
	const Eigen::VectorXd arcIncrement() const;

	// compute arc length of i-th local quintic spline
	double localArc(int i, double t = 1.0, int nqd = 20) const;
	// find internal coordiante t of i-th local spline corresponding given an arclength value
	double arc2t(int i, double arc, double eps = 5.e-14, int nqd = 20) const;
	// i-th order differentiation of spline
	const Eigen::Vector3d d(const Coef &x, int i, double t) const;
	// output spline coefficient to file
	void write(const std::string &name);

private:
	void h(const Eigen::MatrixXd &node);
	void setComponent(int i, Coef &x);
	void computeCoef(Coef &x, BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	static const double qd_GL_x20[20]; // Gauss-Legendre quadrature abscissa
	static const double qd_GL_w20[20]; // Gauss-Legendre quadrature weights
};
} // namespace spline

#endif