#include <eigen3/Eigen/Core>
#include <iostream>

#include "quinticSpline.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI

int main()
{
    // total number of spline segments
    Eigen::Index nodeNum[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    for (Eigen::Index &N : nodeNum)
    {
        std::cout << "Compute global quintic spline of " << std::to_string(N) << " local spline segments ... ";
        // allocate knots input matrix
        Eigen::MatrixX2d knots;
        knots.resize(N + 1, 2);
        // generate knots from analytic shape
        for (int i = 0; i < N + 1; i++)
        {
            double theta = (double)i / N * M_PI;
            knots(i, 0) = sin(theta) * (1 + 0.25 * cos(6 * theta)); // x (or r) coord
            knots(i, 1) = cos(theta) * (1 + 0.25 * cos(6 * theta)); // y (or z) coord
        }
        spline::Quintic sp;
        sp.node(knots);
        // BC of axisymmetric shape
        sp.x(spline::Quintic::BC::Odd, spline::Quintic::BC::Odd);
        sp.y(spline::Quintic::BC::Even, spline::Quintic::BC::Even);

        // write to file
        std::string outputFile = std::to_string(N);
        int zeroNum = 3 - outputFile.length();
        for (int a = 0; a < zeroNum; a++)
        {
            outputFile = "0" + outputFile;
        }
        outputFile = "./resources/output/spline" + outputFile + ".txt";
        std::cout << "write to " << outputFile << std::endl;
        sp.printSpline(outputFile);
    }
    std::cout << "Press enter to finish ... ";
    std::cin.get(); // please stop to let me see
    return 0;
}