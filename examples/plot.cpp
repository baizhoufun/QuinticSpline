#include <eigen3/Eigen/Core>
#include <iostream>

#include "quinticSpline.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

int main()
{
	// total number of spline segments
	Eigen::Index N = 15;
	Eigen::MatrixXd knots;
	knots.resize(N + 1, 2);
	// generate knots from analytic shape
	for (int i = 0; i < N + 1; i++)
	{
		double theta = (double)i / N * M_PI;
		knots(i, 0) = sin(theta) * (1 + 0.26 * cos(4 * theta)); // x
		knots(i, 1) = cos(theta) * (1 + 0.26 * cos(4 * theta)); // y
	}

	spline::Quintic sp;
	sp.init(knots);
	sp.bc(0, spline::Quintic::BCType::Odd, spline::Quintic::BCType::Odd);
	sp.bc(1, spline::Quintic::BCType::Even, spline::Quintic::BCType::Even);

	// BC of axisymmetric shape
	sp.node();
	sp.computeComponent(0);
	sp.computeComponent(1);
	sp.computeArcCoord();

	std::cout << Eigen::MatrixXd(sp.arcCoord()).format(fmt) << "\n";
	int dd = spline::Quintic::search(sp.arcCoord(), sp.arcCoord()[sp.arcCoord().size() - 1]);
	std::cout << dd << "\t";

	return 0;
}