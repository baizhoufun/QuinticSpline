#pragma once
#ifndef SPLINE_H
#define SPLINE_H
#include <eigen3/Eigen/Core>
#include <vector>

namespace spline
{
class Quintic
{

public:
	// enumerate three types of boundary conditions
	enum class BCType
	{
		Even,
		Odd,
		Mix
	};
	struct BC
	{
		BCType bc0, bc1;
		double a0, a1, b0, b1;
	};

private:
	// define data type of quintic spline coefficients, Nx6 matrix
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Coef;
	// ========== data member ==========
	std::vector<Coef> _component;
	Eigen::VectorXd _h;		   // chord
	Eigen::VectorXd _arcCoord; // chord
	Eigen::MatrixXd _node;	   // node_ (input)
	int _dim = 2;
	std::vector<BC> _bc;

public:
	Quintic &operator=(const Quintic &sp);
	Quintic(){};
	Quintic(const Quintic &sp) { *this = sp; }
	void bc(int i, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	const int &dim() const { return _dim; };
	const std::vector<Coef> &component() const { return _component; }; // x-coordinate (or r)
	const Coef &x() const { return _component[0]; };				   // x-coordinate (or r)
	const Coef &y() const { return _component[1]; };				   // x-coordinate (or r)
	const Coef &z() const { return _component[2]; };				   // x-coordinate (or r)
	const Eigen::VectorXd &h() const { return _h; };				   // chord length
	const Eigen::MatrixXd &node() const { return _node; };			   //spline knots
	void init(const Eigen::MatrixXd &xy);
	void node();							 // set node from input matrix
	void node(const Eigen::VectorXd &chord); // set node from input matrix
	void computeComponent(int k);
	// ========== API ==========
	const Eigen::VectorXd &arcCoord() const { return _arcCoord; };

	const Eigen::VectorXd arcIncrement() const;

	// compute arc length of i-th local quintic spline
	double localArc(int i, double t = 1.0, int nqd = 20) const;
	// find internal coordiante t of i-th local spline corresponding given an arclength value
	double arc2t(int i, double arc, double eps = 5.e-14, int nqd = 20) const;
	// i-th order differentiation of spline
	const Eigen::Vector3d d(const Coef &x, int i, double t) const;
	// output spline coefficient to file
	void write(const std::string &name);
	static int search(const Eigen::VectorXd &ar, double key, int low, int high);
	static int search(const Eigen::VectorXd &ar, double key) { return search(ar, key, 0, ar.size() - 1); };
	void computeArcCoord();

private:
	void h(const Eigen::MatrixXd &node);
	void setComponent(int i, Coef &x);
	void computeCoef(Coef &x, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	static const double qd_GL_x20[20]; // Gauss-Legendre quadrature abscissa
	static const double qd_GL_w20[20]; // Gauss-Legendre quadrature weights
};
} // namespace spline

#endif