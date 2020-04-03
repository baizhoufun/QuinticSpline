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
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Coef;
	int _dim = 2;
	std::vector<Coef> _component;
	Eigen::VectorXd _h;		   // chord
	Eigen::VectorXd _arcCoord; // chord
	Eigen::MatrixXd _node;	   // node_ (input)
	std::vector<BC> _bc;

public:
	Quintic &operator=(const Quintic &sp);
	Quintic();
	Quintic(const Quintic &sp);
	const int &dim() const;
	const std::vector<Coef> &component() const; // x-coordinate (or r)
	const Coef &x() const;						// x-coordinate (or r)
	const Coef &y() const;						// x-coordinate (or r)
	const Coef &z() const;						// x-coordinate (or r)
	const Eigen::VectorXd &h() const;			// chord length
	const Eigen::MatrixXd &node() const;		//spline knots
	const Eigen::VectorXd &arcCoord() const;
	// ========== real computation ==========
	void bc(int i, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	void init(const Eigen::MatrixXd &xy);
	void node();							 // set node from input matrix
	void node(const Eigen::VectorXd &chord); // set node from input matrix
	void computeComponent(int k);
	void computeArcCoord();
	// ========== query API ==========
	const Eigen::VectorXd arcIncrement() const;
	double localArc(int i, double t = 1.0, int nqd = 20) const;
	double arc2t(int i, double arc, double eps = 5.e-14, int nqd = 20) const;
	const Eigen::Vector3d d(const Coef &x, int i, double t) const;
	void write(const std::string &name);
	static int search(const Eigen::VectorXd &ar, double key, int low, int high);
	static int search(const Eigen::VectorXd &ar, double key);

private:
	void h(const Eigen::MatrixXd &node);
	void setComponent(int i, Coef &x);
	void computeCoef(Coef &x, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	static const double qd_GL_x20[20]; // Gauss-Legendre quadrature abscissa
	static const double qd_GL_w20[20]; // Gauss-Legendre quadrature weights
};
} // namespace spline

#endif