#include <cmath>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include "quinticSpline.hpp"

namespace spline
{

void Quintic::update()
{
	for (int i = 0; i < dim(); i++)
		computeComponent(i);

	computeArcCoord();
}

int Quintic::search(const Eigen::VectorXd &ar, double key, int low, int high)
{
	//int mid;
	while (low <= high)
	{
		int mid = low + ((high - low) / 2);
		if (ar[mid + 1] >= key && ar[mid] <= key) // key found
			return mid;
		else if (ar[mid] > key) // key may be on the left half
			high = mid - 1;
		else if (ar[mid] < key) // key may be on the right half
			low = mid + 1;
	}
	return -1;
}

Quintic &Quintic::operator=(const Quintic &that)
{
	if (this == &that)
		return *this;
	else
	{
		_h = that.h();
		_arcCoord = that.arcCoord();
		_node = that.node();
		_bc.assign(that.bc().begin(), that.bc().end());
		_component.assign(that.component().begin(), that.component().end());
		_dim = that.dim();
		return *this;
	}
};
void Quintic::setBC(int i, BCType bc0, BCType bc1, double a0, double b0, double a1, double b1)
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

void Quintic::setNode()
{
	h(_node);
	for (int i = 0; i < dim(); i++)
		setComponent(i, _component[i]);
}

void Quintic::setNode(const Eigen::VectorXd &chord)
{
	_h = chord;
	for (int i = 0; i < dim(); i++)
		setComponent(i, _component[i]);
}

void Quintic::h(const Eigen::MatrixXd &node)
{
	_h.setZero(node.rows() - 1);
	for (int i = 0; i < _h.size(); i++)
	{
		double tmp = 0;
		for (int k = 0; k < node.cols(); k++)
			tmp += pow(node(i + 1, k) - node(i, k), 2.0);

		_h(i) = sqrt(tmp);
	}
};

void Quintic::setComponent(int i, CoefficientMat &x)
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
const Eigen::VectorXd Quintic::curvature2D(Curvature2DType curvature2DType) const
{
	Eigen::VectorXd curvature;
	curvature.resize(node().rows());
	curvature.setZero();
	for (Eigen::Index i = 0; i < curvature.size() - 1; i++)
	{
		Eigen::Vector3d xc = d(_component[0], i, 0);
		Eigen::Vector3d yc = d(_component[1], i, 0);

		switch (curvature2DType)
		{
		case Curvature2DType::XY:
		{
			curvature(i) = curvatureXY(xc(1), xc(2), yc(1), yc(2));
			break;
		}
		case Curvature2DType::RZ:
		{
			curvature(i) = curvatureRZ(xc(0), xc(1), xc(2), yc(1), yc(2));
			break;
		}
		default:
			break;
		}
	}
	Eigen::Index last = curvature.size() - 1;
	Eigen::Vector3d xc = d(_component[0], last, 1.0);
	Eigen::Vector3d yc = d(_component[1], last, 1.0);
	switch (curvature2DType)
	{
	case Curvature2DType::XY:
	{
		curvature(last) = curvatureXY(xc(1), xc(2), yc(1), yc(2));
		break;
	}
	case Curvature2DType::RZ:
	{
		curvature(last) = curvatureRZ(xc(0), xc(1), xc(2), yc(1), yc(2));
		break;
	}
	default:
		break;
	}
	return curvature;
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
	/*
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
*/
	computeCoefficient(_component[k], bc0, bc1, a0, b0, a1, b1);
}; // namespace spline

void Quintic::computeCoefficient(CoefficientMat &x, BCType bc0, BCType bc1, double a0, double b0, double a1, double b1)
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

const Eigen::Vector3d Quintic::d(const CoefficientMat &x, int i, double t) const
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
			for (int q = 0; q < dim(); q++)
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

	double f0 = localArc(i, x0, nqd) - arc;
	// find intrinsic coordinate that corresponds to a given arclength fraction - eqn (7.28)
	int counter = 0;
	while (std::abs(f0) > eps)
	{
		double df0 = 0;
		for (int q = 0; q < dim(); q++)
			df0 += pow((d(_component[q], i, x0))(1), 2.0);

		df0 = sqrt(df0);

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

const Eigen::MatrixXd Quintic::arcSample(const Eigen::VectorXd &arcPoints, SampleType sampleType) const
{
	int nA = arcPoints.size();
	Eigen::MatrixXd A(nA, dim());
	A.setZero();
	//double dkey = sp.arcCoord()[sp.arcCoord().size() - 1] / (nA - 1);

	int kBegin = 0, kEnd = 0;

	switch (sampleType)
	{
	case SampleType::HeadSnap:
	{
		A.row(0) = node().row(0);
		kBegin = 1;
		kEnd = nA;
	}
	break;
	case SampleType::TailSnap:
	{
		A.row(nA - 1) = node().row(node().rows() - 1);
		kBegin = 0;
		kEnd = nA - 1;
	}
	break;
	case SampleType::FullSnap:
	{
		A.row(0) = node().row(0);
		A.row(nA - 1) = node().row(node().rows() - 1);
		kBegin = 1;
		kEnd = nA - 1;
	}
	break;
	case SampleType::DontSnap:
	{
		kBegin = 0;
		kEnd = nA;
	}
	break;

	default:
		break;
	}

	for (int k = kBegin; k < kEnd; k++)
	{
		double key = arcPoints(k); //k * dkey * pow((double)k / (nA - 1), 1.1);
		int dd = spline::Quintic::search(arcCoord(), key);
		double t = arc2t(dd, (key - arcCoord()[dd]));
		for (int q = 0; q < dim(); q++)
		{
			A(k, q) = d(component()[q], dd, t)(0);
		}
	}

	return A;
};

void Quintic::write(const std::string &name) const
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

Quintic::Quintic(){};
Quintic::Quintic(const Quintic &that)
{
	_h = that.h();
	_arcCoord = that.arcCoord();
	_node = that.node();
	_bc.assign(that.bc().begin(), that.bc().end());
	_component.assign(that.component().begin(), that.component().end());
	_dim = that.dim();
}
const int &Quintic::dim() const { return _dim; };
const std::vector<Quintic::CoefficientMat> &Quintic::component() const { return _component; }; // x-coordinate (or r)
const Quintic::CoefficientMat &Quintic::x() const { return _component[0]; };				   // x-coordinate (or r)
const Quintic::CoefficientMat &Quintic::y() const
{
	return _component[1 > dim() - 1 ? dim() - 1 : 1];
}; // x-coordinate (or r)
const Quintic::CoefficientMat &Quintic::z() const
{
	return _component[2 > dim() - 1 ? dim() - 1 : 2];
};																// x-coordinate (or r)
const Eigen::VectorXd &Quintic::h() const { return _h; };		// chord length
const Eigen::MatrixXd &Quintic::node() const { return _node; }; //spline knots
const Eigen::VectorXd &Quintic::arcCoord() const { return _arcCoord; };
int Quintic::search(const Eigen::VectorXd &ar, double key) { return search(ar, key, 0, ar.size() - 1); };
const std::vector<Quintic::BC> &Quintic::bc() const { return _bc; }
double Quintic::curvatureXY(double dx, double ddx, double dy, double ddy)
{
	return (dx * ddy - dy * ddx) / pow(dx * dx + dy * dy, 1.5);
}

double Quintic::curvatureRZ(double r, double dr, double ddr, double dz, double ddz)
{
	double eps = 1e-10;

	if (std::abs(r) > eps)
	{
		return (dr * ddz - dz * ddr) / pow(dr * dr + dz * dz, 1.5) + dz / sqrt(dr * dr + dz * dz) / r;
	}
	else
	{
		return 2.0 * ddz / dr / dr;
	}
};

double Quintic::estimateDerivative(int component, BCLocation loc, int derivativeOrder) const
{
	switch (loc)
	{
	case BCLocation::Begin:
	{
		double f0 = node()(0, component);
		double f1 = node()(1, component);
		double f2 = node()(2, component);
		double h0 = h()(0);
		double h1 = h()(1);
		switch (derivativeOrder)
		{
		case 1:
		{
			return (f1 - f0) / h0;
			break;
		}
		case 2:
		{
			return 2.0 * ((f0 - f1) / h0 + (f2 - f1) / h1) / (h0 + h1);
			break;
		}

		default:
			break;
		}
		break;
	}
	case BCLocation::End:
	{
		double f0 = node()(node().rows() - 3, component);
		double f1 = node()(node().rows() - 2, component);
		double f2 = node()(node().rows() - 1, component);
		double h0 = h()(h().size() - 2);
		double h1 = h()(h().size() - 1);
		switch (derivativeOrder)
		{
		case 1:
		{
			return (f2 - f1) / h1;
			break;
		}
		case 2:
		{
			return 2.0 * ((f0 - f1) / h0 + (f2 - f1) / h1) / (h0 + h1);
			break;
		}

		default:
			break;
		}
		break;
	}
	default:
		break;
	}

	return 0;
};
// abscissa and weights of 20-point Gauss-Legendre quadrature rules
const double Quintic::qd_GL_x20[20] = {0.0034357004074525, 0.0180140363610431, 0.0438827858743371, 0.0804415140888906, 0.1268340467699246, 0.1819731596367425, 0.2445664990245865, 0.3131469556422902, 0.3861070744291775, 0.4617367394332513, 0.5382632605667487, 0.6138929255708225, 0.6868530443577098, 0.7554335009754135, 0.8180268403632580, 0.8731659532300750, 0.9195584859111090, 0.9561172141256630, 0.9819859636389570, 0.9965642995925470};
const double Quintic::qd_GL_w20[20] = {0.0088070035695761, 0.0203007149001935, 0.0313360241670545, 0.0416383707883524, 0.0509650599086202, 0.0590972659807592, 0.0658443192245883, 0.0710480546591910, 0.0745864932363019, 0.0763766935653629, 0.0763766935653629, 0.0745864932363019, 0.0710480546591910, 0.0658443192245883, 0.0590972659807592, 0.0509650599086202, 0.0416383707883524, 0.0313360241670545, 0.0203007149001935, 0.0088070035695761};
} // namespace spline