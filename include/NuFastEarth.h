#ifndef NuFastEarth_H
#define NuFastEarth_H

#include <vector>

#include "Earth.h"
#include "Matrix.h"

class Probability_Engine
{
	public:
		Probability_Engine();

		// Oscillation set functions
		void Set_Oscillation_Parameters(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, bool neutrino_mode);
		void Set_s12sq(double s12sq);
		void Set_s13sq(double s13sq);
		void Set_s23sq(double s23sq); // fast
		void Set_delta(double delta); // fast
		void Set_Dmsq21(double Dmsq21);
		void Set_Dmsq31(double Dmsq31);
		void Set_neutrino_mode(bool neutrino_mode);
		// Oscillation get functions
		double Get_s12sq();
		double Get_s13sq();
		double Get_s23sq();
		double Get_delta();
		double Get_Dmsq21();
		double Get_Dmsq31();
		int Get_neutrino_mode_sign();

		void Set_Earth(double detector_depth, Earth_Density *earth_density); // km
		void Set_Production_Height(double production_height); // km, optional, fast
		void Set_rhoYe_Sun(double rhoYe_Sun); // g/cc
		double Get_rhoYe_Sun(); // g/cc
		void Set_Spectra(std::vector<double> Es, std::vector<double> coszs); // GeV

		void Set_Eigenvalue_Precision(int eigenvalue_precision); // not necessary to call as this is set by default

		void Precalc(); // calculates some simple functions that only need to be recalculated when anything other than theta23,delta is changed
		// calculate one set of eigenvalues and eigenvectors (tilde basis)
		Eigen Calculate_Eigen(double rhoYeE);
		Eigen Calculate_Eigen(); // vacuum

		void Calculate_Trajectories();
		void Calculate_Eigens();
		std::vector<std::vector<Eigen>> Get_Eigens_Constant();
		void Calculate_Internal_Amplitudes();
		void Calculate_Eigen_Vac();
		void Calculate_Probabilities();

		// ------------ Solar ------------ //
		void Calculate_Vsolar();
		void Calculate_Solar_Day_In_Earth();
		void Calculate_Solar_Night_In_Earth();

		std::vector<std::vector<Matrix3r>> Get_Probabilities(); // first dimension is energy, second is cosz

		// ------------ Solar ------------ //
		// Does not depend on Earth at all, but does depend on oscillation parameters. The vector is over energies (set in spectra). rhoYe is the density in the Sun in g/cc times the electron fraction
		std::vector<Matrix3r> Get_Solar_Day_Probabilities();
		std::vector<std::vector<Matrix3r>> Get_Solar_Night_Probabilities(); // first dimension is energy, second is cosz. production height should be zero.

	private:
		// input parameters
		double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
		double production_height, detector_depth;
		int neutrino_mode_sign; // +1 for neutrinos, -1 for anti-neutrinos
		Earth_Density *earth_density;
		std::vector<double> Es, coszs;
		double rhoYe_Sun; // density in the Sun, just one point
		int eigenvalue_precision;

		// keep track of whether or not things are set up
		bool oscillation_parameters_set, earth_set, spectra_set, rhoYe_Sun_set;
		// keep track of what has been calculated
		bool trajectories_calculated, oscillation_precalced, eigens_calculated, internal_amplitudes_calculated, eigen_vac_calculated, probabilities_calculated;
		// ------------ Solar ------------ //
		bool Vsolar_calculated, solar_night_in_earth_calculated, solar_day_in_earth_calculated;

		// oscillation physics calculations
		double c12sq, c13sq, s13xc13, A0, B0, Dmsqee, See, Sem, Set, Tee, Tem, Tet;
		double c23sq, s23, c23, s23xc23, s12, c12, c13, s13;
		// calculate one matrix of eigenvectors in the tilde basis
		Matrix3c Probability_Amplitude_1Shell(double LoverE, Eigen evec);
		// does the theta23,delta rotation and computes the probability including adding the production height
		Matrix3r Inner_Amplitude_to_Probability(Matrix3c amp, double cosz, double E);

		// ------------ Solar ------------ //
		// does the theta23,delta rotation and computes the probability from mass state i to flavor state alpha
		Matrix3r Inner_Amplitude_to_Probability_Solar(Matrix3c amp);
		// This computes |V_ai|^2 where V diagonalizes the Hamiltonian in the Sun, takes rhoYe*E*sign
		Matrix3r Calculate_Solar_Probability(double rhoYeE);

		// stored calculations for reuse
		std::vector<std::vector<std::pair<double,double>>> mean_densities1s, mean_densities2s; // holds the trajectories. first index is cosz, second is layers, pair is L,rho
		std::vector<std::vector<Eigen>> eigens_constant; // holds the eigenvalues & eigenvectors (tilde basis) for constant density profiles. first index is energy, second is layers (densities2s first, then densities1s)
		std::vector<std::vector<std::vector<Eigen>>> eigens_varying; // holds the eigenvalues & eigenvectors (tilde basis) for varying density profiles. first index is energy, second is cosz, third is layers (may not all be the same size depending on cosz)
		std::vector<std::vector<Matrix3c>> internal_amplitudes; // holds the internal amplitude matrices (not counting the theta23,delta rotation). first index is energy, second is cosz
		Eigen eigen_vac; // vacuum eigenvalues and eigenvectors in the tilde basis (not counting the theta23,delta rotation)
		std::vector<std::vector<Matrix3r>> probabilities; // first index is energy, second is cosz

		// ------------ Solar ------------ //
		// For solar we won't save the actual probabilities
		std::vector<Matrix3r> Vsolar; // index is energy, flavor to mass basis in the Sun
		Matrix3r solar_day_in_earth; // this is just Usq
		std::vector<std::vector<Matrix3r>> solar_night_in_earth; // first index is energy, second is cosz. mass to flavor basis
};

#endif
