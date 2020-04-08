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
	Eigen::size_t N = 40;
	int dim = 3;
	Eigen::MatrixXd knots;
	knots.resize(N + 1, dim);
	for (int i = 0; i < N + 1; i++)
	{
		double theta = (double)i / N * M_PI;
		knots(i, 0) = 0.9 * theta + 0.6 * sin(4 * theta) * (1 + 0.26 * cos(4 * theta)); // x (or r) coord
		knots(i, 1) = 0.6 * cos(4 * theta) * (1 + 0.26 * cos(4 * theta));				// y (or z) coord
		knots(i, 2) = theta * theta;
	}

	spline::Quintic sp;
	sp.init(knots);
	sp.setNode();

	sp.setBC(0, spline::BCType::Odd, spline::BCType::Odd);
	sp.setBC(1, spline::BCType::Even, spline::BCType::Even);
	double zSlope = (knots(1, 2) - knots(0, 2)) / sp.h()[0];
	sp.setBC(2, spline::BCType::Mix, spline::BCType::Mix, zSlope, 0, zSlope, 0);
	sp.update();

	spline::Quintic sp1(sp);
	for (int g = 0; g < 10; g++)
	{
		sp1.setNode(sp1.arcIncrement());
		sp1.update();
	}
	sp1.write("../examples/01/spline.txt");

	int nA = 30;
	Eigen::VectorXd sample(nA);
	sample.setZero();
	double totalArc = sp1.arcCoord()[sp1.arcCoord().size() - 1];
	for (int k = 0; k < nA; k++)
		sample(k) = totalArc * pow((double)k / (nA - 1), 1.0);

	Eigen::MatrixXd A = sp1.arcSample(sample, spline::SampleType::FullSnap);
	std::ofstream file("../examples/01/scatter.txt");
	file << A.format(fmt) << '\n';
	file.close();

	return 0;
}
