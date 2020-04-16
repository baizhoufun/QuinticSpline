#pragma once
#ifndef QUINTICSPLINE_HPP
#define QUINTICSPLINE_HPP
#include <eigen3/Eigen/Core>
#include <vector>

namespace spline
{
enum class BCType
{
	Even,
	Odd,
	Mix
};
enum class SampleType
{
	HeadSnap,
	TailSnap,
	FullSnap,
	DontSnap
};
enum class Curvature2DType
{
	XY,
	RZ
};

enum class BCLocation
{
	Begin,
	End
};

class Quintic
{
public:
	struct BC
	{
		BCType bc0, bc1;
		double a0, a1, b0, b1;
	};

private:
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> CoefficientMat;
	int _dim = 2;
	std::vector<CoefficientMat> _component;
	Eigen::VectorXd _h;		   // chord
	Eigen::VectorXd _arcCoord; // arc length coordinates
	Eigen::MatrixXd _node;	   // input node
	std::vector<BC> _bc;

public:
	Quintic();								 //constructor
	Quintic(const Quintic &that);			 //copy constructor
	Quintic &operator=(const Quintic &that); //copy assignment
	const int &dim() const;
	const std::vector<CoefficientMat> &component() const; // N-by-6 matrices
	const CoefficientMat &x() const;					  // 1st coordinate
	const CoefficientMat &y() const;					  // 2nd coordinate
	const CoefficientMat &z() const;					  // 3rd coordinate
	const Eigen::VectorXd &h() const;					  // chord length
	const Eigen::MatrixXd &node() const;				  //spline knots
	const Eigen::VectorXd &arcCoord() const;
	const std::vector<BC> &bc() const;
	// ========== real computation ==========
	void init(const Eigen::MatrixXd &xy);
	void setNode();								// set node from input
	void setNode(const Eigen::VectorXd &chord); // set node and chord from input
	void setBC(int i, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	void computeComponent(int k);
	void update();
	// ========== query API ==========
	void write(const std::string &name) const;
	static int search(const Eigen::VectorXd &ar, double key, int low, int high);
	static int search(const Eigen::VectorXd &ar, double key);
	double localArc(int i, double t = 1.0, int nqd = 20) const;
	double arc2t(int i, double arc, double eps = 5.e-14, int nqd = 20) const;
	const Eigen::MatrixXd arcSample(const Eigen::VectorXd &arcPoints, SampleType sampleType = SampleType::DontSnap) const;
	const Eigen::Vector3d d(const CoefficientMat &x, int i, double t) const;
	const Eigen::VectorXd arcIncrement() const;
	const Eigen::VectorXd curvature2D(Curvature2DType curvature2DType) const;
	double estimateDerivative(int component, BCLocation loc, int derivativeOrder) const;

private:
	void h(const Eigen::MatrixXd &node);
	void setComponent(int i, CoefficientMat &x);
	void computeCoefficient(CoefficientMat &x, BCType bc0, BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	void computeArcCoord();

	static const double qd_GL_x20[20]; // Gauss-Legendre quadrature abscissa
	static const double qd_GL_w20[20]; // Gauss-Legendre quadrature weights
	static double curvatureXY(double dx, double ddx, double dy, double ddy);
	static double curvatureRZ(double r, double dr, double ddr, double dz, double ddz);
};
} // namespace spline

#endif