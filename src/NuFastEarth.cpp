#include <vector>
#include <cassert>

#include "NuFastEarth.h"
#include "Earth.h"
#include "Geometry.h"
#include "Matrix.h"

Probability_Engine::Probability_Engine()
{
	oscillation_parameters_set = false;
	earth_set = false;
	spectra_set = false;

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
	trajectories_calculated = false;
	eigens_calculated = false;
	internal_amplitudes_calculated = false;
	probabilities_calculated = false;
	solar_night_in_earth_calculated = false;

	detector_depth = detector_depth_;
	earth_density = earth_density_;

	earth_set = true;
}
void Probability_Engine::Set_Production_Height(double production_height_)
{
	probabilities_calculated = false;
	solar_night_in_earth_calculated = false;

	production_height = production_height_;
}
void Probability_Engine::Set_Spectra(std::vector<double> Es_, std::vector<double> coszs_)
{
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

	// Initialize vectors for the trajectories
	std::vector<std::pair<double,double>> mean_densities1, mean_densities2; // initialize the distances and densities for one cosz path

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
	assert(earth_set);
	assert(spectra_set);

	if (not trajectories_calculated) Calculate_Trajectories();
	if (not oscillation_precalced) Precalc();

	eigens_constant.clear();
	eigens_varying.clear();

	double rhoYe;
	if (earth_density->constant_shells)
	{
		eigens_constant.reserve(Es.size());
		for (unsigned int i = 0; i < Es.size(); i++)
		{
			eigens_constant.emplace_back();
			for (int j = 0; j < earth_density->n_discontinuities; j++)
			{
				eigens_constant[i].emplace_back(Calculate_Eigen(earth_density->rhoYe(earth_density->discontinuities[j] - 1e-8) * Es[i] * neutrino_mode_sign));
			}
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

	eigens_calculated = true;
}
std::vector<std::vector<Eigen>> Probability_Engine::Get_Eigens_Constant()
{
	assert(eigens_calculated);
	assert(earth_density->constant_shells);
	return eigens_constant;
}
void Probability_Engine::Calculate_Internal_Amplitudes()
{
	assert(oscillation_parameters_set);
	assert(earth_set);
	assert(spectra_set);

	if (not trajectories_calculated) Calculate_Trajectories();
	if (not oscillation_precalced) Precalc();
	if (not eigens_calculated) Calculate_Eigens();

	Matrix3c amp;
	Eigen eigen;
	int count;

	internal_amplitudes.clear();

	for (unsigned int i = 0; i < Es.size(); i++)
	{
		internal_amplitudes.emplace_back();
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			if (earth_density->constant_shells)
			{
				amp.Identity();
				count = 0;
				for (unsigned int k = 0; k < mean_densities2s[j].size(); k++)
				{
					eigen = eigens_constant[i][int(mean_densities2s[j][k].second)];
					amp = amp * Probability_Amplitude_1Shell(mean_densities2s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
					count++;
				} // k, mean_densities2s
				amp = amp.AAT(); // double the trajectory
				for (unsigned int k = 0; k < mean_densities1s[j].size(); k++)
				{
					eigen = eigens_constant[i][int(mean_densities1s[j][k].second)];
					amp = amp * Probability_Amplitude_1Shell(mean_densities1s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
					count++;
				}
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
				amp = amp.AAT(); // double the trajectory
				for (unsigned int k = 0; k < mean_densities1s[j].size(); k++)
				{
					eigen = eigens_varying[i][j][count];
					amp = amp * Probability_Amplitude_1Shell(mean_densities1s[j][k].first / Es[i] * neutrino_mode_sign, eigen);
					count++;
				}
				internal_amplitudes[i].emplace_back(amp);
			} // varying shells
		} // j, cosz
	} // i, E

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
	assert(earth_set);
	assert(spectra_set);

	if (not trajectories_calculated) Calculate_Trajectories();
	if (not eigens_calculated) Calculate_Eigens();
	if (not internal_amplitudes_calculated) Calculate_Internal_Amplitudes();
	if (not eigen_vac_calculated) Calculate_Eigen_Vac();

	probabilities.clear();

	probabilities.reserve(Es.size());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		probabilities.emplace_back();
		probabilities[i].reserve(coszs.size());
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			probabilities[i].emplace_back();
			probabilities[i][j] = Inner_Amplitude_to_Probability(internal_amplitudes[i][j], coszs[j], Es[i]);
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
	assert(earth_set);
	assert(spectra_set);

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
