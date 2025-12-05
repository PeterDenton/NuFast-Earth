#include <cmath>
#include <cassert>

#include "NuFastEarth.h"
#include "Matrix.h"
#include "Geometry.h"

#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

constexpr double YerhoE2a = 0.0001526493231029146;
constexpr double eVsqkm_to_GeV_over2 = 1.e-9 / 1.97327e-7 * 1.e3 / 2.; // = 2.5338...

void Probability_Engine::Precalc()
{
	c12sq = 1. - s12sq;
	c13sq = 1. - s13sq;
	s13xc13 = sqrt(s13sq * c13sq);

	A0 = Dmsq21 + Dmsq31;
	B0 = Dmsq21 * Dmsq31;

	Dmsqee = Dmsq31 - s12sq * Dmsq21;

	See =  Dmsq21 + Dmsqee * c13sq;
	Sem = -Dmsq21 * sqrt(c13sq * s12sq * c12sq);
	Set = -Dmsqee * s13xc13;

	Tee =  B0 * c13sq * c12sq;
	Tem =  Dmsq31 * Sem;
	Tet = -B0 * s13xc13 * c12sq;

	c23sq = 1 - s23sq;
	s23 = sqrt(s23sq);
	c23 = sqrt(c23sq);
	s23xc23 = s23 * c23;

	s12 = sqrt(s12sq);
	c12 = sqrt(c12sq);

	s13 = sqrt(s13sq);
	c13 = sqrt(c13sq);

	oscillation_precalced = true;
}

// rhoYeE=rho*Ye*E
// in the tilde basis: does not depend on theta23 or delta
// sign(rhoYeE) indicates nu/nubar
Eigen Probability_Engine::Calculate_Eigen(double rhoYeE)
{
	assert(oscillation_precalced);

	double Amatter, A, B, C, xmat, lambda1, lambda2, lambda3, Dlambda21, rootAsqB, ss0;
	Eigen eigen;

	if (rhoYeE == 0)
	{
		// eigenvalues:
		eigen.lambda.vec[0] = 0;
		eigen.lambda.vec[1] = Dmsq21;
		eigen.lambda.vec[2] = Dmsq31;

		// eigenvectors:
		eigen.Ve2sq = s12sq * c13sq;
		eigen.Ve3sq = s13sq;
		eigen.Vm2sq = c12sq;
		eigen.Vm3sq = 0;
		eigen.Vt2sq = s12sq * s13sq;
		eigen.Vt3sq = c13sq;

		eigen.Ve2Vm2 = s12 * c12 * c13;
		eigen.Ve2Vt2 = -s12sq * c13 * s13;
		eigen.Vm2Vt2 = -c12 * s12 * s13;
		eigen.Ve3Vm3 = 0;
		eigen.Ve3Vt3 = s13 * c13;
		eigen.Vm3Vt3 = 0;

		return eigen;
	} // in vacuum
	else
	{
		// eigenvalues:
		Amatter = YerhoE2a * rhoYeE;

		A = A0 + Amatter;
		B = B0 + Amatter * See;
		C = Amatter * Tee;

		if (eigenvalue_precision >= 0)
		{
			// ---------------------------------- //
			// Get lambda3 from lambda+ of MP/DMP //
			// ---------------------------------- //
			xmat = Amatter / Dmsqee;

			lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + sqrt(sq(1 - xmat) + 4 * s13sq * xmat));
			// ------------------------------------------------------------- //
			// Newton iterations to improve lambda3 arbitrarily, if needed,  //
			// ------------------------------------------------------------- //
			for (int i = 0; i < eigenvalue_precision; i++)
			{
				lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B);
			}
		} // approximate lambda3
		else  // eigenvalue_precision < 0
		{
			// ------------------------------------------- //
			// Get lambda3 from ZS, computationally costly //
			// Lambda3 for both mass orderings			   //
			// ------------------------------------------- //
			rootAsqB = sqrt(sq(A) - 3. * B);
			ss0 = acos((cube(A) - 4.5 * A * B + 13.5 * C) / cube(rootAsqB));
			if (Dmsq31 < 0.)
				ss0 = ss0 + 2. * M_PI; //  add 2Pi if Dmsq31 < 0 to get lambda3  
			lambda3 = (A + 2. * rootAsqB * cos(ss0 / 3.)) / 3.;
		} // calculate lambda3 exactly

		// ---------------------------------------------------- //
		// Use A,C to get other eigenvalues and  Delta lambda's //
		// ---------------------------------------------------- //
		Dlambda21 = sqrt(sq(A - lambda3) - 4. * C / lambda3);
		if (Dmsq21 < 0.) Dlambda21 = -Dlambda21; // relevant for checking solar mass ordering
		lambda2 = 0.5 * (A - lambda3 + Dlambda21);
		lambda1 = lambda2 - Dlambda21;

		eigen.lambda.vec[0] = lambda1;
		eigen.lambda.vec[1] = lambda2;
		eigen.lambda.vec[2] = lambda3;

		// eigenvectors:
		double Dlambda21, Dlambda31;
		double InvPID, Xp2, Xp3, InvV3xV3, InvV2xV2;

		Dlambda21 = lambda2 - lambda1;
		Dlambda31 = lambda3 - lambda1;

		InvPID = 1. / (Dlambda31 * (Dlambda31 - Dlambda21) * Dlambda21);

		Xp2 = -Dlambda31 * InvPID;
		Xp3 =  Dlambda21 * InvPID;

		eigen.Ve3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
		eigen.Ve2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

		eigen.Ve3Vm3 = (-lambda3 * Sem + Tem) * Xp3;
		eigen.Ve3Vt3 = (-lambda3 * Set + Tet) * Xp3;

		InvV3xV3 = eigen.Ve3Vm3 / eigen.Ve3sq;
		eigen.Vm3Vt3 = eigen.Ve3Vt3 * InvV3xV3;

		eigen.Vm3sq = eigen.Ve3Vm3 * InvV3xV3;
		eigen.Vt3sq = 1. - eigen.Ve3sq - eigen.Vm3sq;

		eigen.Ve2Vm2 = (-lambda2 * Sem + Tem) * Xp2;
		eigen.Ve2Vt2 = (-lambda2 * Set + Tet) * Xp2;

		InvV2xV2 = eigen.Ve2Vm2 / eigen.Ve2sq;
		eigen.Vm2Vt2 = eigen.Ve2Vt2 * InvV2xV2;

		eigen.Vm2sq = eigen.Ve2Vm2 * InvV2xV2;
		eigen.Vt2sq = 1. - eigen.Ve2sq - eigen.Vm2sq;

		return eigen;
	} // in matter
}
Eigen Probability_Engine::Calculate_Eigen() // vacuum
{
	return Calculate_Eigen(0);
}

// LoverE = L/E
Matrix3c Probability_Engine::Probability_Amplitude_1Shell(double LoverE, Eigen evec)
{
	double Dlambda31, Dlambda21;
	double D31, D21, cosD31, sinD31, cosD21, sinD21, cosD31m1, cosD21m1;
	Matrix3c prob_amp;

	Dlambda31 = evec.lambda.vec[2] - evec.lambda.vec[0];
	Dlambda21 = evec.lambda.vec[1] - evec.lambda.vec[0];

	D31 =  eVsqkm_to_GeV_over2 * Dlambda31 * LoverE;
	D21 =  eVsqkm_to_GeV_over2 * Dlambda21 * LoverE;

	cosD31 = cos(D31);
	sinD31 = sin(D31);

	cosD21 = cos(D21);
	sinD21 = sin(D21);

	cosD31m1 = cosD31 - 1.;
	cosD21m1 = cosD21 - 1.;

	// do the diagonal first
	prob_amp.arr[0][0] = std::complex<double>(1 + evec.Ve2sq * cosD21m1 + evec.Ve3sq * cosD31m1, -evec.Ve2sq * sinD21 - evec.Ve3sq * sinD31);
	prob_amp.arr[1][1] = std::complex<double>(1 + evec.Vm2sq * cosD21m1 + evec.Vm3sq * cosD31m1, -evec.Vm2sq * sinD21 - evec.Vm3sq * sinD31);
	prob_amp.arr[2][2] = std::complex<double>(1 + evec.Vt2sq * cosD21m1 + evec.Vt3sq * cosD31m1, -evec.Vt2sq * sinD21 - evec.Vt3sq * sinD31);

	prob_amp.arr[0][1] = std::complex<double>(evec.Ve2Vm2 * cosD21m1 + evec.Ve3Vm3 * cosD31m1, -evec.Ve2Vm2 * sinD21 - evec.Ve3Vm3 * sinD31);
	prob_amp.arr[0][2] = std::complex<double>(evec.Ve2Vt2 * cosD21m1 + evec.Ve3Vt3 * cosD31m1, -evec.Ve2Vt2 * sinD21 - evec.Ve3Vt3 * sinD31);
	prob_amp.arr[1][2] = std::complex<double>(evec.Vm2Vt2 * cosD21m1 + evec.Vm3Vt3 * cosD31m1, -evec.Vm2Vt2 * sinD21 - evec.Vm3Vt3 * sinD31);

	prob_amp.arr[1][0] = prob_amp.arr[0][1];
	prob_amp.arr[2][0] = prob_amp.arr[0][2];
	prob_amp.arr[2][1] = prob_amp.arr[1][2];

	return prob_amp;
}

// does the theta23,delta rotation and computes the probability
Matrix3r Probability_Engine::Inner_Amplitude_to_Probability(Matrix3c amp, double cosz, double E)
{
	Matrix3r prob;
	double sinzsq, L, xtmp, ytmp;
	std::complex<double> epid, emid, Ate, Aet, AmtpAtm;

	if (production_height > 0)
	{
		sinzsq = 1 - sq(cosz);
		L = Distance_Between_Layers(earth_density->r_E + production_height, earth_density->r_E, sinzsq, detector_depth, earth_density);
		amp = amp * Probability_Amplitude_1Shell(L / E, eigen_vac);
	}

	epid = exp(std::complex<double>(0, delta)); // e^(+i*delta)
	emid = conj(epid); // e^(-i*delta)

	Ate = epid * amp.arr[2][0];
	Aet = emid * amp.arr[0][2];

	AmtpAtm = emid * amp.arr[1][2] +  epid * amp.arr[2][1];

	xtmp = amp.arr[0][0].real();
	ytmp = amp.arr[0][0].imag();
	prob.arr[0][0] = sq(xtmp) + sq(ytmp);

	xtmp = c23 * amp.arr[1][0].real() + s23 * Ate.real();
	ytmp = c23 * amp.arr[1][0].imag() + s23 * Ate.imag();
	prob.arr[0][1] = sq(xtmp) + sq(ytmp);

	prob.arr[0][2] = 1. - prob.arr[0][0] - prob.arr[0][1];

	xtmp = c23 * amp.arr[0][1].real() + s23 * Aet.real();
	ytmp = c23 * amp.arr[0][1].imag() + s23 * Aet.imag();
	prob.arr[1][0] = sq(xtmp) + sq(ytmp);

	xtmp = c23sq * amp.arr[1][1].real() + s23sq * amp.arr[2][2].real() + s23xc23 * AmtpAtm.real();
	ytmp = c23sq * amp.arr[1][1].imag() + s23sq * amp.arr[2][2].imag() + s23xc23 * AmtpAtm.imag();
	prob.arr[1][1] = sq(xtmp) + sq(ytmp);

	prob.arr[1][2]= 1. - prob.arr[1][0] - prob.arr[1][1];

	prob.arr[2][0]= 1. - prob.arr[0][0] - prob.arr[1][0];
	prob.arr[2][1]= 1. - prob.arr[0][1] - prob.arr[1][1];
	prob.arr[2][2]= 1. - prob.arr[0][2] - prob.arr[1][2];

	return prob;
}

// Takes an inner amplitude (with no theta23 and delta rotation) and converts it to a matrix of probabilities for solar neutrinos
// Instead of alpha to beta, it is mass state i to beta
Matrix3r Probability_Engine::Inner_Amplitude_to_Probability_Solar(Matrix3c amp)
{
	Matrix3r prob;

	std::complex<double> eid;
	eid = exp(std::complex<double>(0, delta));

	// fill in the probability matrix
	prob.arr[0][0] = std::norm(amp.arr[0][0] * c12 * c13 - amp.arr[0][1] * s12 - amp.arr[0][2] * c12 * s13);
	prob.arr[0][1] = std::norm(amp.arr[0][0] * s12 * c13 + amp.arr[0][1] * c12 - amp.arr[0][2] * s12 * s13);
	prob.arr[0][2] = 1 - prob.arr[0][0] - prob.arr[0][1];

	prob.arr[1][1] = std::norm(c12 * (amp.arr[1][1] * c23 + amp.arr[2][1] * eid * s23) + s12 * (c13 * (amp.arr[1][0] * c23 + amp.arr[2][0] * eid * s23) - s13 * (amp.arr[1][2] * c23 + amp.arr[2][2] * eid * s23)));
	prob.arr[1][2] = std::norm(s13 * (amp.arr[1][0] * c23 + amp.arr[2][0] * eid * s23) + c13 * (amp.arr[1][2] * c23 + amp.arr[2][2] * eid * s23));
	prob.arr[1][0] = 1 - prob.arr[1][1] - prob.arr[1][2];

	prob.arr[2][0] = 1 - prob.arr[0][0] - prob.arr[1][0];
	prob.arr[2][1] = 1 - prob.arr[0][1] - prob.arr[1][1];
	prob.arr[2][2] = 1 - prob.arr[0][2] - prob.arr[1][2];

	return prob.Transpose();
}
// This computes |V_ai|^2 where V diagonalizes the Hamiltonian in the Sun
// takes rhoYe*E*sign
Matrix3r Probability_Engine::Calculate_Solar_Probability(double rhoYeE)
{
	// Get the eigens
	Eigen eigen = Calculate_Eigen(rhoYeE);

	double cosdelta;
	Matrix3r prob;

	cosdelta = cos(delta);

	// fill in the array in the real flavor basis
	prob.arr[0][1] = eigen.Ve2sq;
	prob.arr[0][2] = eigen.Ve3sq;
	prob.arr[0][0] = 1 - prob.arr[0][1] - prob.arr[0][2];

	prob.arr[1][1] = c23sq * eigen.Vm2sq + s23sq * eigen.Vt2sq + 2 * cosdelta * s23xc23 * eigen.Vm2Vt2;
	prob.arr[1][2] = c23sq * eigen.Vm3sq + s23sq * eigen.Vt3sq + 2 * cosdelta * s23xc23 * eigen.Vm3Vt3;
	prob.arr[1][0] = 1 - prob.arr[1][1] - prob.arr[1][2];

	prob.arr[2][0] = 1 - prob.arr[1][0] - prob.arr[0][0];
	prob.arr[2][1] = 1 - prob.arr[1][1] - prob.arr[0][1];
	prob.arr[2][2] = 1 - prob.arr[1][2] - prob.arr[0][2];

	return prob;
}
