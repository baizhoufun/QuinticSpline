#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <iostream>

#include "quinticSpline.hpp"

namespace spline
{
int Quintic::search(const Eigen::VectorXd &ar, double key, int low, int high)
{
	int mid;
	while (low < high - 1)
	{
		mid = low + ((high - low) / 2);
		if (ar[mid + 1] >= key && ar[mid] <= key)
			return mid;
		if (ar[mid] > key) // key may be on the left half
			high = mid;
		else if (ar[mid] < key) // key may be on the right half
			low = mid;
	}
	return -1;
}

void Quintic::write(const std::string &name)
{
	Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
	std::ofstream file(name);
	for (size_t i = 0; i < component().size(); i++)
	{
		file << component()[i].format(fmt) << '\n';
	}
	file << h().format(fmt) << '\n';
	file.close();
}

Quintic &Quintic::operator=(const Quintic &sp)
{
	if (this == &sp)
		return *this;
	else
	{
		// 	_node.resize(sp.node().rows(), 2);
		// 	_node = sp.node();
		// 	_h.resize(sp.h().size());
		// 	_h = sp.h();
		// 	_x.resize(sp.x().rows(), 6);
		// 	_x = sp.x();
		// 	_y.resize(sp.y().rows(), 6);
		// 	_y = sp.y();
		// 	return *this;
	}
};
void Quintic::bc(int i, BCType bc0, BCType bc1, double a0, double b0, double a1, double b1)
{
	_bc[i].bc0 = bc0;
	_bc[i].bc1 = bc1;
	_bc[i].a0 = a0;
	_bc[i].a1 = a1;
	_bc[i].b0 = b0;
	_bc[i].b1 = b1;
};

void Quintic::init(const Eigen::MatrixXd &xy)
{
	_node.resize(xy.rows(), xy.cols());
	_node = xy;
	_dim = xy.cols();
	_component.resize(dim());
	_bc.resize(dim());
}

void Quintic::node()
{
	h(_node);
	for (Eigen::size_t i = 0; i < dim(); i++)
	{
		setComponent(i, _component[i]);
	}
	//computeArcCoord();
}

void Quintic::node(const Eigen::VectorXd &chord)
{
	_h = chord;
	for (Eigen::size_t i = 0; i < dim(); i++)
	{
		setComponent(i, _component[i]);
	}
	//computeArcCoord();
}

void Quintic::h(const Eigen::MatrixXd &node)
{
	_h.setZero(node.rows() - 1);
	for (Eigen::size_t i = 0; i < _h.size(); i++)
	{
		double tmp = 0;
		for (Eigen::size_t k = 0; k < node.cols(); k++)
			tmp += pow(node(i + 1, k) - node(i, k), 2.0);

		_h(i) = sqrt(tmp);
	}
};

void Quintic::setComponent(int i, Coef &x)
{
	x.resize(_node.rows(), 6);
	x.setZero();
	x.col(0) = _node.col(i);
};
void Quintic::computeArcCoord()
{

	_arcCoord.resize(h().size() + 1);
	_arcCoord.setZero();
	for (Eigen::Index i = 1; i < _arcCoord.size(); i++)
	{
		_arcCoord[i] = _arcCoord[i - 1] + localArc(i - 1);
	}
};
const Eigen::VectorXd Quintic::arcIncrement() const
{
	Eigen::VectorXd tmp;
	tmp.resize(h().size());
	tmp.setZero();
	for (Eigen::Index i = 0; i < tmp.size(); i++)
	{
		tmp[i] = localArc(i);
	}
	return tmp;
};

void Quintic::computeComponent(int k)
{

	const BCType &bc0 = _bc[k].bc0;
	const BCType &bc1 = _bc[k].bc1;

	const double &a0 = _bc[k].a0;
	const double &a1 = _bc[k].a1;
	const double &b0 = _bc[k].b0;
	const double &b1 = _bc[k].b1;

	printf("comp. %d BC = ", k);
	switch (bc0)
	{
	case BCType::Odd:
		printf("odd ");
		break;
	case BCType::Even:
		printf("even ");
		break;
	case BCType::Mix:
		printf("mix ");
		break;
	default:
		break;
	}

	switch (bc1)
	{
	case BCType::Odd:
		printf("odd\n");
		break;
	case BCType::Even:
		printf("even\n");
		break;
	case BCType::Mix:
		printf("mix\n");
		break;
	default:
		break;
	}

	computeCoef(_component[k], bc0, bc1, a0, b0, a1, b1);
}; // namespace spline

void Quintic::computeCoef(Coef &x, BCType bc0, BCType bc1, double a0, double b0, double a1, double b1)
{
	Eigen::Index N = _node.rows() - 1;
	Eigen::VectorXd rhs(2 * N + 2);
	const Eigen::VectorXd X = x.col(0);
	rhs.setZero();
	Eigen::SparseMatrix<double> A(2 * N + 2, 2 * N + 2);
	A.reserve(Eigen::VectorXi::Constant(A.cols(), 6));

	switch (bc0)
	{
	case BCType::Even:
	{ // even condition - eqn (7.17)
		rhs(0) = 0;
		rhs(1) = -10. * X(0) + 10. * X(1);
		double h0 = _h(0);
		A.insert(0, 0) = 1.0;
		A.insert(1, 0) = 6. * h0;
		A.insert(1, 1) = 3. * h0 * h0;
		A.insert(1, 2) = 4. * h0;
		A.insert(1, 3) = -h0 * h0;
		break;
	}
	case BCType::Odd:
	{ // even condition - eqn (7.18)
		rhs(0) = 15. * X(0) - 15. * X(1);
		rhs(1) = 0;
		double h0 = _h(0);
		A.insert(0, 0) = -8. * h0;
		A.insert(0, 1) = -3. * h0 * h0;
		A.insert(0, 2) = -7. * h0;
		A.insert(0, 3) = 2. * h0 * h0;
		A.insert(1, 1) = 1.0;
		break;
	}
	case BCType::Mix:
	{ // mix condition - eqn (7.19)
		rhs(0) = a0;
		rhs(1) = b0 / 2;
		A.insert(0, 0) = 1.0;
		A.insert(1, 1) = 1.0;
		break;
	}
	default:
		break;
	}
	switch (bc1)
	{
	case BCType::Even:
	{ // even condition - eqn (7.23)
		rhs(rhs.size() - 2) = 0;
		rhs(rhs.size() - 1) = 10. * X(X.size() - 2) - 10. * X(X.size() - 1);
		double hm1 = _h(_h.size() - 1);
		A.insert(A.rows() - 2, A.rows() - 2) = 1.0;

		A.insert(A.rows() - 1, A.rows() - 4) = -4. * hm1;
		A.insert(A.rows() - 1, A.rows() - 3) = -1. * hm1 * hm1;
		A.insert(A.rows() - 1, A.rows() - 2) = -6. * hm1;
		A.insert(A.rows() - 1, A.rows() - 1) = +3. * hm1 * hm1;
		break;
	}
	case BCType::Odd:
	{ // even condition - eqn (7.24)
		rhs(rhs.size() - 2) = 15. * X(X.size() - 2) - 15. * X(X.size() - 1);
		rhs(rhs.size() - 1) = 0;
		double hm1 = _h(_h.size() - 1);
		A.insert(A.rows() - 2, A.rows() - 4) = -7. * hm1;
		A.insert(A.rows() - 2, A.rows() - 3) = -2. * hm1 * hm1;
		A.insert(A.rows() - 2, A.rows() - 2) = -8. * hm1;
		A.insert(A.rows() - 2, A.rows() - 1) = +3. * hm1 * hm1;
		A.insert(A.rows() - 1, A.rows() - 1) = 1.0;
		break;
	}
	case BCType::Mix:
	{ // mix condition - eqn (7.25)
		rhs(rhs.size() - 2) = a1;
		rhs(rhs.size() - 1) = b1 / 2.;
		A.insert(A.rows() - 2, A.rows() - 2) = 1.0;
		A.insert(A.rows() - 1, A.rows() - 1) = 1.0;
		break;
	}
	default:
		break;
	}

	// spline coefficient c1 and c2 matrix entries - eqn (7.13)
	for (Eigen::Index j = 1; j <= N - 1; j++)
	{
		double lb = _h(j) / _h(j - 1);
		double lb2 = lb * lb;
		double lb3 = lb * lb * lb;
		double lb4 = lb * lb * lb * lb;
		double hj = _h(j);
		double hj2 = _h(j) * _h(j);

		rhs(2 * j) = -10. * (lb3 * X(j - 1) - (1. + lb3) * X(j) + X(j + 1));
		rhs(2 * j + 1) = 15. * (-lb4 * X(j - 1) + (lb4 - 1.) * X(j) + X(j + 1));

		A.insert(2 * j, 2 * j - 2) = 4.0 * hj * lb2;
		A.insert(2 * j, 2 * j - 1) = 1.0 * hj2 * lb;
		A.insert(2 * j, 2 * j) = 6.0 * hj * (lb2 - 1.);
		A.insert(2 * j, 2 * j + 1) = -3.0 * hj2 * (1. + lb);
		A.insert(2 * j, 2 * j + 2) = -4.0 * hj;
		A.insert(2 * j, 2 * j + 3) = 1.0 * hj2;

		A.insert(2 * j + 1, 2 * j - 2) = 7.0 * hj * lb3;
		A.insert(2 * j + 1, 2 * j - 1) = 2.0 * hj2 * lb2;
		A.insert(2 * j + 1, 2 * j) = 8.0 * hj * (1. + lb3);
		A.insert(2 * j + 1, 2 * j + 1) = 3.0 * hj2 * (1. - lb2);
		A.insert(2 * j + 1, 2 * j + 2) = 7.0 * hj;
		A.insert(2 * j + 1, 2 * j + 3) = -2.0 * hj2;
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
	solver.compute(A);
	Eigen::VectorXd lhs = solver.solve(rhs);

	Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2>> c1(lhs.data(), lhs.size() / 2);
	Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2>> c2(lhs.data() + 1, lhs.size() / 2);

	x.col(1) = c1;
	x.col(2) = c2;

	// reconstruction of c3 c4 c5 - eqn (7.10 - 7.12)
	for (Eigen::Index j = 0; j <= N - 1; j++)
	{
		double hj = _h(j);
		double c1j = c1(j), c2j = c2(j), c1jp1 = c1(j + 1), c2jp1 = c2(j + 1);
		// re-parametrize local spline - eqn (7.27)
		x(j, 1) = c1j * hj;
		x(j, 2) = c2j * hj * hj;
		x(j, 3) = -6. * hj * c1j - 3. * hj * hj * c2j - 4. * hj * c1jp1 + 1. * hj * hj * c2jp1 + 10. * (X(j + 1) - X(j));
		x(j, 4) = +8. * hj * c1j + 3. * hj * hj * c2j + 7. * hj * c1jp1 - 2. * hj * hj * c2jp1 - 15. * (X(j + 1) - X(j));
		x(j, 5) = -3. * hj * c1j - 1. * hj * hj * c2j - 3. * hj * c1jp1 + 1. * hj * hj * c2jp1 + 6. * (X(j + 1) - X(j));
	}
};

const Eigen::Vector3d Quintic::d(const Coef &x, int i, double t) const
{
	double t2 = t * t;
	double t3 = t * t * t;
	double t4 = t * t * t * t;
	double t5 = t * t * t * t * t;
	double x0 = x(i, 0);
	double x1 = x(i, 1);
	double x2 = x(i, 2);
	double x3 = x(i, 3);
	double x4 = x(i, 4);
	double x5 = x(i, 5);
	Eigen::Vector3d tmp(0, 0, 0);
	// derivatives of local quintic spline with respect to t
	tmp(0) = x0 + t * x1 + t2 * x2 + t3 * x3 + t4 * x4 + t5 * x5;
	tmp(1) = x1 + 2. * t * x2 + 3. * t2 * x3 + 4. * t3 * x4 + 5. * t4 * x5;
	tmp(2) = 2. * x2 + 6. * t * x3 + 12. * t2 * x4 + 20. * t3 * x5;
	return tmp;
};

double Quintic::localArc(int i, double t, int nqd) const
{
	// local arclength of each local spline - eqn (7.28)
	if (i < _node.rows() - 1 && i >= 0)
	{
		// numerical integration using 20-point Gauss-Legendre quadrature rules
		const double *qdx = qd_GL_x20;
		const double *qdw = qd_GL_w20;
		double arc = 0;
		for (int k = 0; k < nqd; k++)
		{
			double ab = t * qdx[k];
			double tmp = 0;
			for (Eigen::size_t q = 0; q < dim(); q++)
				tmp += pow((d(component()[q], i, ab))(1), 2.0);

			arc += qdw[k] * sqrt(tmp);
		}
		return t * arc;
	}
	else
	{
		printf("Spline index out of range! ");
		return 0.;
	}
};

double Quintic::arc2t(int i, double arc, double eps, int nqd) const
{
	double x0 = 0.5;
	double epsilon = eps;
	double f0 = localArc(i, x0, nqd) - arc;
	// find intrinsic coordiante that corresponds to a given arclength fraction - eqn (7.28)
	int counter = 0;
	while (std::abs(f0) > epsilon)
	{
		double df0 = 0;
		for (Eigen::size_t q = 0; q < dim(); q++)
			df0 += pow((d(_component[q], i, x0))(1), 2.0);

		//= sqrt(pow((d(_x, i, x0))(1), 2.0) + pow((d(_y, i, x0))(1), 2.0));
		x0 = x0 - f0 / df0;
		f0 = localArc(i, x0, nqd) - arc;
		counter++;
		if (counter > 10)
		{
			printf("arc2t not converging .. try lowering error tolerance");
		}
	}
	return x0;
};

// abscissa and weights of 20-point Gauss-Legendre quadrature rules
const double Quintic::qd_GL_x20[20] = {0.0034357004074525, 0.0180140363610431, 0.0438827858743371, 0.0804415140888906, 0.1268340467699246, 0.1819731596367425, 0.2445664990245865, 0.3131469556422902, 0.3861070744291775, 0.4617367394332513, 0.5382632605667487, 0.6138929255708225, 0.6868530443577098, 0.7554335009754135, 0.8180268403632580, 0.8731659532300750, 0.9195584859111090, 0.9561172141256630, 0.9819859636389570, 0.9965642995925470};
const double Quintic::qd_GL_w20[20] = {0.0088070035695761, 0.0203007149001935, 0.0313360241670545, 0.0416383707883524, 0.0509650599086202, 0.0590972659807592, 0.0658443192245883, 0.0710480546591910, 0.0745864932363019, 0.0763766935653629, 0.0763766935653629, 0.0745864932363019, 0.0710480546591910, 0.0658443192245883, 0.0590972659807592, 0.0509650599086202, 0.0416383707883524, 0.0313360241670545, 0.0203007149001935, 0.0088070035695761};
} // namespace spline