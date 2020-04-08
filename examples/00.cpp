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
        Eigen::MatrixX2d knots;
        knots.resize(N + 1, 2);

        for (int i = 0; i < N + 1; i++)
        {
            double theta = (double)i / N * M_PI;
            knots(i, 0) = sin(theta) * (1 + 0.25 * cos(8 * theta)); // x (or r) coord
            knots(i, 1) = cos(theta) * (1 + 0.25 * cos(8 * theta)); // y (or z) coord
        }

        spline::Quintic sp;
        sp.init(knots);
        sp.setNode(); // BC of axis-symmetric shape
        sp.setBC(0, spline::BCType::Odd, spline::BCType::Odd);
        sp.setBC(1, spline::BCType::Even, spline::BCType::Even);
        sp.update();

        // write to file
        std::string outputFile = std::to_string(N);
        int zeroNum = 3 - outputFile.length();
        for (int a = 0; a < zeroNum; a++)
            outputFile = "0" + outputFile;

        outputFile = "../examples/00/spline" + outputFile + ".txt";
        std::cout << "write to " << outputFile << std::endl;
        sp.write(outputFile);
    }

    return 0;
}