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
	Eigen::Index nodeNum[] = {20};
	Eigen::Index N = 200;
	Eigen::MatrixXd knots;
	knots.resize(N + 1, 2);
	// generate knots from analytic shape
	for (int i = 0; i < N + 1; i++)
	{
		double theta = (double)i / N * M_PI;
		knots(i, 0) = sin(theta) * (1 + 0.26 * cos(4 * theta)); // x (or r) coord
		knots(i, 1) = cos(theta) * (1 + 0.26 * cos(4 * theta)); // y (or z) coord
	}
	spline::Quintic sp;
	sp.node(knots);

	// BC of axisymmetric shape
	sp.computeComponent(0, spline::Quintic::BC::Odd, spline::Quintic::BC::Odd);
	sp.computeComponent(1, spline::Quintic::BC::Even, spline::Quintic::BC::Even);
	std::cout << "=========== h ===========\n";
	std::cout << Eigen::MatrixXd(sp.h()).format(fmt) << "\n";
	std::cout << "=========== arc ===========\n";
	std::cout << Eigen::MatrixXd(sp.arcIncrement()).format(fmt) << "\n";
	std::cout << "=========== arcCoord ===========\n";
	std::cout << Eigen::MatrixXd(sp.arcCoord()).format(fmt) << "\n";

	for (int g = 0; g < 30; g++)
	{

		std::cout << "Compute global quintic spline of " << std::to_string(g) << " local spline segments ... \n";
		// allocate knots input matrix

		//spline::Quintic sp1;
		sp.node(knots, sp.arcIncrement());

		// BC of axisymmetric shape
		sp.computeComponent(0, spline::Quintic::BC::Odd, spline::Quintic::BC::Odd);
		sp.computeComponent(1, spline::Quintic::BC::Even, spline::Quintic::BC::Even);

		// write to file
		// std::string outputFile = std::to_string(N);
		// int zeroNum = 3 - outputFile.length();
		// for (int a = 0; a < zeroNum; a++)
		// {
		// 	outputFile = "0" + outputFile;
		// }
		// outputFile = "./resources/output/spline" + outputFile + ".txt";
		// std::cout << "write to " << outputFile << std::endl;
		// sp.write(outputFile);

		std::cout << "=========== h ===========\n";
		std::cout << Eigen::MatrixXd(sp.h()).format(fmt) << "\n";
		std::cout << "=========== arc ===========\n";
		std::cout << Eigen::MatrixXd(sp.arcIncrement()).format(fmt) << "\n";
		std::cout << "=========== arcCoord ===========\n";
		std::cout << Eigen::MatrixXd(sp.arcCoord()).format(fmt) << "\n";
	}

	std::cout << "Press enter to finish ... ";
	std::cin.get(); // please stop to let me see
	return 0;
}