#include <eigen3/Eigen/Core>
#include <iostream>
#include <fstream>

#include "quinticSpline.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

int main()
{
	// total number of spline segments
	Eigen::Index N = 50;
	int dim = 3;
	Eigen::MatrixXd knots;
	knots.resize(N + 1, dim);
	// generate knots from analytic shape
	for (int i = 0; i < N + 1; i++)
	{
		double theta = (double)i / N * M_PI;
		knots(i, 0) = 0.9 * theta + 0.6 * sin(4 * theta) * (1 + 0.26 * cos(4 * theta)); // x (or r) coord
		knots(i, 1) = 0.6 * cos(4 * theta) * (1 + 0.26 * cos(4 * theta));				// y (or z) coord
		knots(i, 2) = theta;
	}

	spline::Quintic sp;
	sp.init(knots);
	sp.node();

	sp.bc(0, spline::BCType::Odd, spline::BCType::Odd);
	sp.bc(1, spline::BCType::Even, spline::BCType::Even);
	double zslope = (knots(1, 2) - knots(0, 2)) / sp.h()[0];
	sp.bc(2, spline::BCType::Mix, spline::BCType::Mix, zslope, 0, zslope, 0);

	sp.update();

	for (int g = 0; g < 10; g++)
	{
		sp.node(sp.arcIncrement());
		sp.update();
	}

	sp.write("../resources/spline.txt");

	int nA = 10;
	Eigen::MatrixXd A(nA, dim);
	A.setZero();
	double dkey = sp.arcCoord()[sp.arcCoord().size() - 1] / (nA - 1);

	for (int k = 0; k < nA; k++)
	{
		double key = k * dkey;
		int dd = spline::Quintic::search(sp.arcCoord(), key);
		double t = sp.arc2t(dd, (key - sp.arcCoord()[dd]));
		A(k, 0) = sp.d(sp.x(), dd, t)(0);
		A(k, 1) = sp.d(sp.y(), dd, t)(0);
		A(k, 2) = sp.d(sp.z(), dd, t)(0);
	}

	std::ofstream file("../resources/scatter.txt");
	file << A.format(fmt) << '\n';
	file.close();

	return 0;
}
