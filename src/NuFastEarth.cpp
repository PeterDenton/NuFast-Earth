#include <vector>
#include <cassert>

#include "NuFastEarth.h"
#include "Earth.h"
#include "Geometry.h"
#include "Matrix.h"

namespace NuFast {

Probability_Engine::Probability_Engine()
{
	oscillation_parameters_set = false;
	earth_set = false;
	spectra_set = false;
	E_spectra_set = false;
	trajectory_set = false;

	earth_mode = false;
	single_trajectory_mode = false;

	trajectories_calculated = false;
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false; // the ones not counting the theta23,delta rotation
	eigen_vac_calculated = false;

	probabilities_calculated = false;

	production_height = 0; // default to zero km

	// ------------ Solar ------------ //
	rhoYe_Sun_set = false;

	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	eigenvalue_precision = 1; // default to 1 Newton-Raphson iteration, set to negative for exact eigenvalues
}

void Probability_Engine::Set_Oscillation_Parameters(double s12sq_, double s13sq_, double s23sq_, double delta_, double Dmsq21_, double Dmsq31_, bool neutrino_mode_)
{
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	eigen_vac_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	s12sq = s12sq_;
	s13sq = s13sq_;
	s23sq = s23sq_;
	delta = delta_;
	Dmsq21 = Dmsq21_;
	Dmsq31 = Dmsq31_;
	neutrino_mode_sign = neutrino_mode_ ? +1 : -1;

	oscillation_parameters_set = true;
}

void Probability_Engine::Set_s12sq(double s12sq_)
{
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	eigen_vac_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	s12sq = s12sq_;
}

void Probability_Engine::Set_s13sq(double s13sq_)
{
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	eigen_vac_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	s13sq = s13sq_;
}

void Probability_Engine::Set_s23sq(double s23sq_)
{
	oscillation_precalced = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	s23sq = s23sq_;
}

void Probability_Engine::Set_delta(double delta_)
{
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;
	solar_day_in_earth_calculated = false;

	delta = delta_;
}

void Probability_Engine::Set_Dmsq21(double Dmsq21_)
{
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	eigen_vac_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	Dmsq21 = Dmsq21_;
}

void Probability_Engine::Set_Dmsq31(double Dmsq31_)
{
	oscillation_precalced = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	eigen_vac_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	Dmsq31 = Dmsq31_;
}

void Probability_Engine::Set_neutrino_mode(bool neutrino_mode_)
{
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	neutrino_mode_sign = neutrino_mode_ ? +1 : -1;
}

double Probability_Engine::Get_s12sq()
{
	assert(oscillation_parameters_set);
	return s12sq;
}
double Probability_Engine::Get_s13sq()
{
	assert(oscillation_parameters_set);
	return s13sq;
}
double Probability_Engine::Get_s23sq()
{
	assert(oscillation_parameters_set);
	return s23sq;
}
double Probability_Engine::Get_Dmsq21()
{
	assert(oscillation_parameters_set);
	return Dmsq21;
}
double Probability_Engine::Get_Dmsq31()
{
	assert(oscillation_parameters_set);
	return Dmsq31;
}
double Probability_Engine::Get_delta()
{
	assert(oscillation_parameters_set);
	return delta;
}
int Probability_Engine::Get_neutrino_mode_sign()
{
	assert(oscillation_parameters_set);
	return neutrino_mode_sign;
}
void Probability_Engine::Set_Earth(double detector_depth_, Earth_Density *earth_density_)
{
	assert(not single_trajectory_mode);

	trajectories_calculated = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	solar_night_in_earth_calculated = false;

	detector_depth = detector_depth_;
	earth_density = earth_density_;

	earth_set = true;
	earth_mode = true;
}
void Probability_Engine::Set_Production_Height(double production_height_)
{
	assert(not single_trajectory_mode);

	probabilities_calculated = false;
	solar_night_in_earth_calculated = false;

	production_height = production_height_;
}
void Probability_Engine::Set_Spectra(std::vector<double> Es_, std::vector<double> coszs_)
{
	assert(not single_trajectory_mode);

	trajectories_calculated = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	// check that every element is positive
	for (unsigned int i = 0; i < Es_.size(); i++)
		assert(Es_[i] > 0);

	Es = Es_;
	coszs = coszs_;

	spectra_set = true;
	earth_mode = true;
}
void Probability_Engine::Set_E_Spectra(std::vector<double> Es_)
{
	assert(not earth_mode);

	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	// check that every element is positive
	for (unsigned int i = 0; i < Es_.size(); i++)
		assert(Es_[i] > 0);

	Es = Es_;

	E_spectra_set = true;
	single_trajectory_mode = true;
}
void Probability_Engine::Set_Trajectory(std::vector<std::pair<double,double>> trajectory_)
{
	assert(not earth_mode);

	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	// check the elements
	for (unsigned int i = 0; i < trajectory_.size(); i++)
	{
		assert(trajectory_[i].first > 0); // length, km
		assert(trajectory_[i].second >= 0); // rhoYe, g/cc
	} // i, trajectory

	trajectory = trajectory_;

	trajectory_set = true;
	single_trajectory_mode = true;
}
void Probability_Engine::Set_rhoYe_Sun(double rhoYe_Sun_) // g/cc
{
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	rhoYe_Sun = rhoYe_Sun_;

	rhoYe_Sun_set = true;
}
double Probability_Engine::Get_rhoYe_Sun()
{
	assert(rhoYe_Sun_set);
	return rhoYe_Sun;
}
void Probability_Engine::Set_Eigenvalue_Precision(int eigenvalue_precision_)
{
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	Vsolar_calculated = false;
	solar_night_in_earth_calculated = false;

	eigenvalue_precision = eigenvalue_precision_;
}

void Probability_Engine::Calculate_Trajectories()
{
	assert(spectra_set);
	assert(earth_mode);

	// Initialize vectors for the trajectories
	std::vector<std::pair<double,double>> mean_densities1, mean_densities2; // initialize the distances and densities for one cosz path

	mean_densities1s.clear();
	mean_densities2s.clear();

	mean_densities1s.reserve(coszs.size());
	mean_densities2s.reserve(coszs.size());

	// Get the trajectories
	for (unsigned int i = 0; i < coszs.size(); i++)
	{
		Mean_Densities(coszs[i], detector_depth, earth_density, mean_densities1, mean_densities2);
		mean_densities1s.emplace_back(mean_densities1);
		mean_densities2s.emplace_back(mean_densities2);
	} // i, looping through the cosz vector

	trajectories_calculated = true;
}

void Probability_Engine::Calculate_Eigens()
{
	assert(oscillation_parameters_set);

	if (earth_mode)
	{
		assert(earth_set);
		assert(spectra_set);
		if (not trajectories_calculated) Calculate_Trajectories();
	}
	else if (single_trajectory_mode)
	{
		assert(E_spectra_set);
		assert(trajectory_set);
	}
	else
		assert(false);

	if (not oscillation_precalced) Precalc();

	eigens_constant.clear();
	eigens_varying.clear();

	double rhoYe;

	if (earth_mode)
	{
		if (earth_density->constant_shells)
		{
			eigens_constant.reserve(Es.size());
			for (unsigned int i = 0; i < Es.size(); i++)
			{
				eigens_constant.emplace_back();
				for (int j = 0; j < earth_density->n_discontinuities; j++)
					eigens_constant[i].emplace_back(Calculate_Eigen(earth_density->rhoYe(earth_density->discontinuities[j] - 1e-8) * Es[i] * neutrino_mode_sign));
			} // i, E
		} // constant shells
		else
		{
			for (unsigned int i = 0; i < Es.size(); i++)
			{
				eigens_varying.emplace_back();
				for (unsigned int j = 0; j < coszs.size(); j++)
				{
					eigens_varying[i].emplace_back();
					for (unsigned int k = 0; k < mean_densities2s[j].size(); k++)
					{
						rhoYe = mean_densities2s[j][k].second;
						eigens_varying[i][j].emplace_back(Calculate_Eigen(rhoYe * Es[i] * neutrino_mode_sign));
					} // k, trajectory2
					for (unsigned int k = 0; k < mean_densities1s[j].size(); k++) // loop over layers
					{
						rhoYe = mean_densities1s[j][k].second;
						eigens_varying[i][j].emplace_back(Calculate_Eigen(rhoYe * Es[i] * neutrino_mode_sign));
					} // k, trajectory1
				} // j, coszs
			} // i, E
		} // varying shells
	} // earth_mode
	else if (single_trajectory_mode)
	{
		eigens_constant.reserve(Es.size());
		for (unsigned int i = 0; i < Es.size(); i++)
		{
			eigens_constant.emplace_back();
			for (unsigned int j = 0; j < trajectory.size(); j++)
				eigens_constant[i].emplace_back(Calculate_Eigen(trajectory[j].second * Es[i] * neutrino_mode_sign));
		} // i, E
	}

	eigens_calculated = true;
}
std::vector<std::vector<Eigen>> Probability_Engine::Get_Eigens_Constant()
{
	assert(eigens_calculated);
	if (earth_mode)
		assert(earth_density->constant_shells);
	return eigens_constant;
}
void Probability_Engine::Calculate_Internal_Amplitudes()
{
	assert(oscillation_parameters_set);
	if (earth_mode)
	{
		assert(earth_set);
		assert(spectra_set);
		if (not trajectories_calculated) Calculate_Trajectories();
	}
	else if (single_trajectory_mode)
	{
		assert(trajectory_set);
		assert(E_spectra_set);
	}
	else
		assert(false);

	if (not oscillation_precalced) Precalc();
	if (not eigens_calculated) Calculate_Eigens();

	Matrix3c amp;
	Eigen eigen;
	int count;

	internal_amplitudes.clear();

	if (earth_mode)
	{
		for (unsigned int i = 0; i < Es.size(); i++)
		{
			internal_amplitudes.emplace_back();
			for (unsigned int j = 0; j < coszs.size(); j++)
			{
				if (earth_density->constant_shells)
				{
					amp.Identity();
					for (unsigned int k = 0; k < mean_densities2s[j].size(); k++)
					{
						eigen = eigens_constant[i][int(mean_densities2s[j][k].second)];
						amp = amp * Probability_Amplitude_1Shell(mean_densities2s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
					} // k, mean_densities2s
					// We have now calculated A_O
					amp = amp.AAT(); // double the trajectory
					// We have now calculated A_O*A_I
					for (int k = mean_densities1s[j].size() - 1; k >=0; k--)
					{
						eigen = eigens_constant[i][int(mean_densities1s[j][k].second)];
						amp = amp * Probability_Amplitude_1Shell(mean_densities1s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
					} // k, mean_densities1s
					// We have now calculated A_O*A_I*A_S
					internal_amplitudes[i].emplace_back(amp);
				} // constant shells
				else
				{
					amp.Identity();
					count = 0;
					for (unsigned int k = 0; k < mean_densities2s[j].size(); k++)
					{
						eigen = eigens_varying[i][j][count];
						amp = amp * Probability_Amplitude_1Shell(mean_densities2s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
						count++;
					} // k, mean_densities2s
					// We have now calculated A_O
					amp = amp.AAT(); // double the trajectory
					// We have now calculated A_O*A_I
					for (int k = mean_densities1s[j].size() - 1; k >= 0; k--)
					{
						eigen = eigens_varying[i][j][count];
						amp = amp * Probability_Amplitude_1Shell(mean_densities1s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
						count++;
					}
					// We have now calculated A_O*A_I*A_S
					internal_amplitudes[i].emplace_back(amp);
				} // varying shells
			} // j, cosz
		} // i, E
	} // earth_mode
	else if (single_trajectory_mode)
	{
		for (unsigned int i = 0; i < Es.size(); i++)
		{
			internal_amplitudes.emplace_back();
			amp.Identity();
			for (unsigned int j = 0; j < trajectory.size(); j++)
			{
				eigen = eigens_constant[i][j];
				amp = Probability_Amplitude_1Shell(trajectory[j].first / Es[i] * neutrino_mode_sign, eigen) * amp;
			} // j, trajectory
			internal_amplitudes[i].emplace_back(amp);
		} // i, E
	} // single trajectory mode

	internal_amplitudes_calculated = true;
}
void Probability_Engine::Calculate_Eigen_Vac()
{
	assert(oscillation_parameters_set);

	eigen_vac = Calculate_Eigen();

	eigen_vac_calculated = true;
}
void Probability_Engine::Calculate_Probabilities()
{
	assert(oscillation_parameters_set);
	int n_cosz;
	if (earth_mode)
	{
		assert(earth_set);
		assert(spectra_set);
		if (not trajectories_calculated) Calculate_Trajectories();
		n_cosz = coszs.size();
	}
	else if (single_trajectory_mode)
	{
		assert(E_spectra_set);
		assert(trajectory_set);
		n_cosz = 1;
	}
	else
		assert(false);

	if (not eigens_calculated) Calculate_Eigens();
	if (not internal_amplitudes_calculated) Calculate_Internal_Amplitudes();
	if (not eigen_vac_calculated) Calculate_Eigen_Vac();

	probabilities.clear();

	probabilities.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		probabilities.emplace_back();
		probabilities[i].reserve(n_cosz);
		for (int j = 0; j < n_cosz; j++)
		{
			probabilities[i].emplace_back();
			if (earth_mode)
				probabilities[i][j] = Inner_Amplitude_to_Probability(internal_amplitudes[i][j], coszs[j], Es[i]);
			else if (single_trajectory_mode)
				probabilities[i][j] = Inner_Amplitude_to_Probability(internal_amplitudes[i][j], 0, Es[i]);
		} // j, coszs
	} // i, Es

	probabilities_calculated = true;
}

// ------------ Solar ------------ //
void Probability_Engine::Calculate_Vsolar()
{
	assert(oscillation_parameters_set);
	assert(spectra_set);

	if (not oscillation_precalced) Precalc();

	Vsolar.clear();

	Vsolar.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		Vsolar.emplace_back();
		Vsolar[i] = Calculate_Solar_Probability(rhoYe_Sun * Es[i] * neutrino_mode_sign);
	} // i, Es

	Vsolar_calculated = true;
}
void Probability_Engine::Calculate_Solar_Night_In_Earth()
{
	assert(oscillation_parameters_set);
	assert(earth_set);
	assert(spectra_set);

	if (not oscillation_precalced) Precalc();
	if (not trajectories_calculated) Calculate_Trajectories();
	if (not eigens_calculated) Calculate_Eigens();
	if (not internal_amplitudes_calculated) Calculate_Internal_Amplitudes();

	solar_night_in_earth.clear();

	solar_night_in_earth.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		solar_night_in_earth.emplace_back();
		solar_night_in_earth[i].reserve(coszs.size());
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			solar_night_in_earth[i].emplace_back();
			solar_night_in_earth[i][j] = Inner_Amplitude_to_Probability_Solar(internal_amplitudes[i][j]);
		} // j, coszs
	} // i, Es

	solar_night_in_earth_calculated = true;
}
void Probability_Engine::Calculate_Solar_Day_In_Earth()
{
	assert(oscillation_parameters_set);
	assert(spectra_set);
	assert(rhoYe_Sun_set);

	if (not oscillation_precalced) Precalc();

	// precalc some stuff
	double cosdelta = cos(delta);

	// calculate the norm square of elements of the pmns matrix
	Matrix3r Usq, UsqT;
	Usq.arr[0][0] = c12sq * c13sq;
	Usq.arr[0][1] = s12sq * c13sq;
	Usq.arr[0][2] = s13sq;

	Usq.arr[1][0] = s12sq * c23sq + c12sq * s23sq * s13sq + 2 * cosdelta * s12 * c12 * s23 * c23 * s13;
	Usq.arr[1][2] = s23sq * c13sq;
	Usq.arr[1][1] = 1 - Usq.arr[1][0] - Usq.arr[1][2];

	Usq.arr[2][0] = 1 - Usq.arr[0][0] - Usq.arr[1][0];
	Usq.arr[2][1] = 1 - Usq.arr[0][1] - Usq.arr[1][1];
	Usq.arr[2][2] = 1 - Usq.arr[0][2] - Usq.arr[1][2];

	UsqT = Usq.Transpose();

	solar_day_in_earth = UsqT;

	solar_day_in_earth_calculated = true;
}
std::vector<std::vector<Matrix3r>> Probability_Engine::Get_Probabilities()
{
	assert(oscillation_parameters_set);
	if (earth_mode)
	{
		assert(earth_set);
		assert(spectra_set);
	}
	else if (single_trajectory_mode)
	{
		assert(E_spectra_set);
		assert(trajectory_set);
	}
	else
		assert(false);

	if (not probabilities_calculated) Calculate_Probabilities();
	return probabilities;
}
// first dimension is energy
std::vector<Matrix3r> Probability_Engine::Get_Solar_Day_Probabilities()
{
	assert(oscillation_parameters_set);
	assert(rhoYe_Sun_set);
	assert(spectra_set);
	if (not Vsolar_calculated) Calculate_Vsolar();
	if (not solar_day_in_earth_calculated) Calculate_Solar_Day_In_Earth();

	std::vector<Matrix3r> probabilities_solar_day;

	probabilities_solar_day.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		probabilities_solar_day.emplace_back();
		probabilities_solar_day[i] = Vsolar[i] * solar_day_in_earth;
	} // i, Es

	return probabilities_solar_day;
}
// first dimension is energy, second is cosz. production height should be zero.
std::vector<std::vector<Matrix3r>> Probability_Engine::Get_Solar_Night_Probabilities()
{
	assert(oscillation_parameters_set);
	assert(earth_set);
	assert(production_height == 0);
	assert(rhoYe_Sun_set);
	assert(spectra_set);

	if (not Vsolar_calculated) Calculate_Vsolar();
	if (not solar_night_in_earth_calculated) Calculate_Solar_Night_In_Earth();

	std::vector<std::vector<Matrix3r>> probabilities_solar_night;

	probabilities_solar_night.clear();

	probabilities_solar_night.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		probabilities_solar_night.emplace_back();
		probabilities_solar_night[i].reserve(coszs.size());
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			probabilities_solar_night[i].emplace_back();
			probabilities_solar_night[i][j] = Vsolar[i] * solar_night_in_earth[i][j];
		} // j, coszs
	} // i, Es

	return probabilities_solar_night;
}

// ------------ Print, I/O ------------ //
void Probability_Engine::Print_Status()
{
	printf("\nProbability engine current status:\n");

	printf("\tProduction height = %g km\n", production_height);
	printf("\tEigenvalue precision = %d (0: good, 1: great, ..., negative=exact)\n", eigenvalue_precision);

	if (oscillation_parameters_set)
	{
		printf("\tOscillation parameters have been set.\n");
		printf("\t\ts12sq = %g\n", s12sq);
		printf("\t\ts13sq = %g\n", s13sq);
		printf("\t\ts23sq = %g\n", s23sq);
		printf("\t\tdelta = %g rad\n", delta);
		printf("\t\tDmsq21 = %g eV^2\n", Dmsq21);
		printf("\t\tDmsq31 = %g eV^2\n", Dmsq31);
	} // oscillation parameters set
	else
		printf("\tOscillation parameters have not been set.\n");

	if (earth_mode)
	{
		printf("\tEarth mode selected.\n");

		if (earth_set)
		{
			printf("\tEarth model has been set.\n");
			printf("\t\tDetector depth = %g km\n", detector_depth);
		} // earth set
		else
			printf("\tEarth model has not been set.\n");

		if (spectra_set)
		{
			printf("\tSpectra have been set.\n");
			printf("\t\tNumber of energy grid points = %lu\n", Es.size());
			printf("\t\tNumber of cosz grid points = %lu\n", coszs.size());

		} // spectra set
		else
			printf("\tSpectra have not been set.\n");

	} // earth mode
	else if (single_trajectory_mode)
	{
		printf("\tSingle trajectory mode selected.\n");

		if (trajectory_set)
		{
			printf("\tTrajectory has been set.\n");
			printf("\t\tTrajectory has %lu steps.\n", trajectory.size());
		} // trajectory set
		else
			printf("\tTrajectory has not been set.\n");

		if (E_spectra_set)
		{
			printf("\tEnergy spectrum has been set.\n");
			printf("\t\tNumber of energy grid points = %lu\n", Es.size());
		} // E spectra set
		else
			printf("\tEnergy spectrum has not been set.\n");
	} // single trajectory mode
	else
		printf("\tNo trajectory mode selected.\n");

	if (rhoYe_Sun_set)
		printf("\tDensity in the Sun has been set.\n");
	else
		printf("\tDensity in the Sun has not been set.\n");

	if (earth_mode)
	{
		if (trajectories_calculated)
			printf("\tTrajectories have been calculated.\n");
		else
			printf("\tTrajectories have not been calculated.\n");
	}

	if (oscillation_precalced)
		printf("\tOscillation parameter precalculations have been calculated.\n");
	else
		printf("\tOscillation parameter precalculations have not been calculated.\n");

	if (eigens_calculated)
		printf("\tEigenvalue and eigenvector information has been calculated.\n");
	else
		printf("\tEigenvalue and eigenvector information has not been calculated.\n");

	if (eigen_vac_calculated)
		printf("\tVacuum eigenvalue and eigenvector information has been calculated.\n");
	else
		printf("\tVacuum eigenvalue and eigenvector information has not been calculated.\n");

	if (internal_amplitudes_calculated)
		printf("\tInternal amplitudes have been calculated.\n");
	else
		printf("\tInternal amplitudes have not been calculated.\n");

	if (probabilities_calculated)
		printf("\tProbabilities have been calculated.\n");
	else
		printf("\tProbabilities have not been calculated.\n");
	
	// ------------ Solar ------------ //
	if (Vsolar_calculated)
		printf("\tSolar eigenvector information has been calculated.\n");
	else
		printf("\tSolar eigenvector information has not been calculated.\n");

	if (solar_day_in_earth_calculated)
		printf("\tSolar day time detection information has been calculated.\n");
	else
		printf("\tSolar day time detection information has not been calculated.\n");

	if (solar_night_in_earth_calculated)
		printf("\tSolar night time inner amplitude has been calculated.\n");
	else
		printf("\tSolar night time inner amplitude has not been calculated.\n");

	printf("\n");
}

} // namespace NuFast
