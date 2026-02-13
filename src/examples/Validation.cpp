#include <cassert>

#include "examples/Validation.h"
#include "NuFastEarth.h"
#include "Matrix.h"
#include "Geometry.h"
#include "Solar.h"

void Eigenvalues_Validation()
{
	Probability_Engine probability_engine;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6
	probability_engine.Set_Eigenvalue_Precision(-1);
	double rho, Ye;
	rho = 3;
	Ye = 0.5;
	Constant earth_density(rho * Ye);
	probability_engine.Set_Earth(0, &earth_density);

	double E_min, E_max, E_step;
	int n;
	n = 1e3;
	E_min = 1e-5;
	E_max = 15;
	E_step = (E_max - E_min) / n;
	std::vector<double> Es;
	for (int i = 0; i <= n; i++)
		Es.emplace_back(E_min + i * E_step);
	std::vector coszs = {-0.1};
	probability_engine.Set_Spectra(Es, coszs);

	FILE *data = fopen("data/Eigenvalues_Validation.txt", "w");
	fprintf(data, "%g %g\n", rho, Ye);
	probability_engine.Set_neutrino_mode(false);
	probability_engine.Calculate_Eigens();
	std::vector<std::vector<Eigen>> eigens_constant = probability_engine.Get_Eigens_Constant();
	for (int i = (int)Es.size() - 1; i > -1; i--)
	{
		fprintf(data, "%g ", -Es[i]);
		for (int j = 0; j < 3; j++)
			fprintf(data, "%g ", eigens_constant[i][0].lambda.vec[j]);
		fprintf(data, "\n");
	} // i, E
	probability_engine.Set_neutrino_mode(true);
	probability_engine.Calculate_Eigens();
	eigens_constant = probability_engine.Get_Eigens_Constant();
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		fprintf(data, "%g ", Es[i]);
		for (int j = 0; j < 3; j++)
			fprintf(data, "%g ", eigens_constant[i][0].lambda.vec[j]);
		fprintf(data, "\n");
	} // i, E
	fclose(data);
}

void LBL_Appearance_Validation(bool normal_ordering, bool lower_octant)
{
	// Initialize probability engine
	Probability_Engine probability_engine;

	// Set initial oscillation parameters
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6

	// Set mass ordering
	if (not normal_ordering)	probability_engine.Set_Dmsq31(-2.534e-3);
	if (lower_octant)			probability_engine.Set_s23sq(0.439);

	// Create Earth model instance
	Constant earth_density(3 * 0.5);

	// Set Earth details
	probability_engine.Set_Earth(0, &earth_density);

	// Initialize the energy and cosz arrays
	double E_min, E_max, E_step;
	int n;
	n = 1e3;
	E_min = 1.0;
	E_max = 5;
	E_step = (E_max - E_min) / n;
	std::vector<double> Es;
	for (int i = 0; i <= n; i++)
		Es.emplace_back(E_min + i * E_step);
	std::vector coszs = {L2cosz(1297)}; // DUNE baseline
	probability_engine.Set_Spectra(Es, coszs);

	// Compute the probabilities for different values of delta
	std::vector<std::vector<std::vector<Matrix3r>>> prob_nus, prob_nubars;
	for (int i = 0; i < 4; i++)
	{
		probability_engine.Set_delta(i * M_PI / 2);
		probability_engine.Set_neutrino_mode(true);
		prob_nus.emplace_back(probability_engine.Get_Probabilities());
		probability_engine.Set_neutrino_mode(false);
		prob_nubars.emplace_back(probability_engine.Get_Probabilities());
	}

	std::string fname = "data/LBL_Appearance_Validation_";
	fname += normal_ordering ? "NO" : "IO";
	fname += "_";
	fname += lower_octant ? "lower" : "upper";
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");

	for (unsigned int i = 0; i < coszs.size(); i++)
		fprintf(data, "%g ", coszs[i]);
	fprintf(data, "\n");
	fprintf(data, "%g %g\n", probability_engine.Get_s23sq(), probability_engine.Get_Dmsq31());

	for (unsigned int i = 0; i < Es.size(); i++)
	{
		fprintf(data, "%g ", Es[i]);
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			for (int k = 0; k < 4; k++)
			{
				fprintf(data, "%g ", prob_nus[k][i][j].arr[1][0]);
			} // k, 4, delta
			for (int k = 0; k < 4; k++)
			{
				fprintf(data, "%g ", prob_nubars[k][i][j].arr[1][0]);
			} // k, 4, delta
		} // j, Es
		fprintf(data, "\n");
	} // i, cosz
	fclose(data);
}
void LBL_Appearance_Validation()
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			LBL_Appearance_Validation(i == 0, j == 0);
	} // i, ordering, 2
}
void Solar_Day_Validation()
{
	// Initialize probability engine
	Probability_Engine probability_engine;

	// Set initial oscillation parameters
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6

	// Create Earth model instance
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);

	// Set Earth details
	probability_engine.Set_Earth(2, &earth_density);

	probability_engine.Set_rhoYe_Sun(100 * (2. / 3));

	double E_min, E_max, E_scale;
	int n;

	n = 1e3;

	E_min = 0.1e-3;
	E_max = 50e-3;
	E_scale = pow(E_max / E_min, 1. / n);

	std::vector<double> Es;
	for (int i = 0; i <= n; i++)
		Es.emplace_back(E_min * pow(E_scale, i));

	std::vector<double> coszs = {1.};
	probability_engine.Set_Spectra(Es, coszs);

	// Compute the probabilities
	std::vector<Matrix3r> probs;
	probs = probability_engine.Get_Solar_Day_Probabilities();

	FILE *data = fopen("data/Solar_Day_Validation.txt", "w");
	fprintf(data, "%g %g %g %g\n", probability_engine.Get_Dmsq21(), probability_engine.Get_s12sq(), probability_engine.Get_s13sq(), probability_engine.Get_rhoYe_Sun());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		fprintf(data, "%g %g\n", Es[i], probs[i].arr[0][0]);
	} // i, E
	fclose(data);
}
void Solar_Day_Night_Validation()
{
	// Initialize probability engine
	Probability_Engine probability_engine;

	// Set initial oscillation parameters
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6

	// Create Earth model instance
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);

	// Set Earth details
	// Detector depth
	probability_engine.Set_Earth(2, &earth_density);

	probability_engine.Set_rhoYe_Sun(100 * (2. / 3));

	std::vector<double> Es, coszs;
	double E_min, E_max, E_step, cosz_min, cosz_max, cosz_step;
	int n;

	n = 1e3;

	E_min = 5e-3;
	E_max = 20e-3;
	E_step = (E_max - E_min) / n;

	cosz_min = -1;
	cosz_max = -0.0;
	cosz_step = (cosz_max - cosz_min) / n;

	Es.reserve(n + 1);
	coszs.reserve(n + 1);
	for (int i = 0; i <= n; i++)
	{
		Es.emplace_back(E_min + i * E_step);
		coszs.emplace_back(cosz_min + i * cosz_step);
	} // i, n

	probability_engine.Set_Spectra(Es, coszs);

	// Compute the probabilities
	std::vector<std::vector<Matrix3r>> night_probs;
	std::vector<Matrix3r> day_probs;
	night_probs = probability_engine.Get_Solar_Night_Probabilities();
	day_probs = probability_engine.Get_Solar_Day_Probabilities();

	FILE *data = fopen("data/Solar_Day_Night_Validation.txt", "w");
	double eta;
	double latitude_DUNE, weight_DUNE, sum_DUNE, sum_weight_DUNE;
	double latitude_SK, weight_SK, sum_SK, sum_weight_SK;
	latitude_DUNE = 44.356 * M_PI / 180; // DUNE far detector latitude, see end of table 1 in 1707.02322
	latitude_SK = 36.48 * M_PI / 180; // kamioka, see fig. 3 hep-ph/9702343
	fprintf(data, "%g %g %g %g\n", probability_engine.Get_Dmsq21(), probability_engine.Get_s12sq(), probability_engine.Get_s13sq(), probability_engine.Get_rhoYe_Sun());
	for (unsigned int i = 0; i < Es.size(); i++)
	{
		sum_DUNE = 0;
		sum_SK = 0;
		sum_weight_DUNE = 0;
		sum_weight_SK = 0;
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			eta = M_PI - acos(coszs[j]) + 1e-8;
			weight_DUNE = Solar_Weight(eta, latitude_DUNE) / fabs(sin(eta));
			weight_SK = Solar_Weight(eta, latitude_SK) / fabs(sin(eta));
			sum_DUNE += night_probs[i][j].arr[0][0] * weight_DUNE;
			sum_SK += night_probs[i][j].arr[0][0] * weight_SK;
			sum_weight_DUNE += weight_DUNE;
			sum_weight_SK += weight_SK;
		}
		fprintf(data, "%g %g %g %g\n", Es[i], day_probs[i].arr[0][0], sum_DUNE / sum_weight_DUNE, sum_SK / sum_weight_SK);
	} // i, E
	fclose(data);
}
void Solar_Weight_Validation()
{
	double latitude_DUNE, latitude_SK;
	double eta, eta_min, eta_max, eta_step;
	int n;

	latitude_DUNE = 44.356 * M_PI / 180; // DUNE far detector latitude, see end of table 1 in 1707.02322
	latitude_SK = 36.48 * M_PI / 180; // kamioka, see fig. 3 hep-ph/9702343

	n = 1e3;
	eta_min = 0;
	eta_max = 0.5 * M_PI;
	eta_step = (eta_max - eta_min) / n;

	FILE *data = fopen("data/Solar_Weight_Validation.txt", "w");
	for (int i = 0; i <= n; i++)
	{
		eta = eta_min + i * eta_step;
		fprintf(data, "%g %g %g\n", eta, Solar_Weight(eta, latitude_DUNE), Solar_Weight(eta, latitude_SK));
	} // i, eta, n
	fclose(data);
}
void Diff(std::vector<std::vector<Matrix3r>> probs1, std::vector<std::vector<Matrix3r>> probs2, int alpha1, int beta1, int alpha2, int beta2, double *mean_diff, double *std_diff, double *max_diff)
{
	assert(probs1.size() == probs2.size());
	assert(probs1[0].size() == probs2[0].size());
	assert(alpha1 >= 0 and alpha1 < 3);
	assert(beta1 >= 0 and beta1 < 3);
	assert(alpha2 >= 0 and alpha2 < 3);
	assert(beta2 >= 0 and beta2 < 3);

	double Pab1_diff, Pab2_diff, sum_diff, sum_diff_sq;
	int n;

	sum_diff = 0;
	sum_diff_sq = 0;
	*max_diff = -1;
	for (unsigned int i = 0; i < probs1.size(); i++)
	{
		for (unsigned int j = 0; j < probs1[0].size(); j++)
		{
			Pab1_diff = fabs(probs1[i][j].arr[alpha1][beta1] - probs2[i][j].arr[alpha1][beta1]);
			Pab2_diff = fabs(probs1[i][j].arr[alpha2][beta2] - probs2[i][j].arr[alpha2][beta2]);

			sum_diff += Pab1_diff;
			sum_diff += Pab2_diff;

			sum_diff_sq += sq(Pab1_diff);
			sum_diff_sq += sq(Pab2_diff);

			*max_diff = fmax(*max_diff, Pab1_diff);
			*max_diff = fmax(*max_diff, Pab2_diff);
		} // j
	} // i

	n = int(2 * probs1.size() * probs1[0].size());
	*mean_diff = sum_diff / n;
	*std_diff = sqrt(sum_diff_sq / n - sq(*mean_diff));
}
std::vector<std::vector<Matrix3r>> Precision(int eigenvalue_precision, int n_layer)
{
	Probability_Engine probability_engine;

	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6
	probability_engine.Set_Eigenvalue_Precision(eigenvalue_precision);

	Earth_Density *earth_density;
	if (abs(n_layer + 1) < 1e-3)
		earth_density = new PREM_Full();
	else if (abs(n_layer + 2) < 1e-3)
		earth_density = new PREM_Four();
	else
	{
		int n_discontinuity_layer = int(0.25 * n_layer);
		if (4 * n_discontinuity_layer < n_layer) n_discontinuity_layer++;
		earth_density = new PREM_NDiscontinuityLayer(n_discontinuity_layer);
	}

	probability_engine.Set_Earth(2, earth_density);
	probability_engine.Set_Production_Height(10);

	// Create the energy and cosz vectors
	std::vector<double> Es, coszs;
	double Emin, Emax, Estep, coszmin, coszmax, coszstep;
	int n;

	n = 1e2;

	Emin = 2;
	Emax = 40;

	coszmin = -1;
	coszmax = +0.1;

	Estep = (Emax - Emin) / n;
	coszstep = (coszmax - coszmin) / n;

	Es.reserve(n + 1);
	coszs.reserve(n + 1);
	for (int i = 0; i <= n; i++)
	{
		Es.emplace_back(Emin + i * Estep);
		coszs.emplace_back(coszmin + i * coszstep);
	} // i, n

	probability_engine.Set_Spectra(Es, coszs);

	std::vector<std::vector<Matrix3r>> probabilities = probability_engine.Get_Probabilities();
	delete earth_density;
	return probabilities;

}
void Precision()
{
	int eigenvalue_precision, n_layer, n_layer_true, n_discontinuity_layer;
	double mean_diff, std_diff, max_diff;
	std::vector<std::vector<Matrix3r>> probs_true, probs_test;

	n_layer_true = 1e6; // how many layers to get to the ``true'' probabilities

	probs_true = Precision(-1, n_layer_true);

	FILE *data = fopen("data/Precision.txt", "w");

	for (eigenvalue_precision = 0; eigenvalue_precision <= 2; eigenvalue_precision++)
	{
		printf("Eigenvalue precision: %i\n", eigenvalue_precision);
		probs_test = Precision(eigenvalue_precision, n_layer_true);
		Diff(probs_true, probs_test, 1, 1, 1, 0, &mean_diff, &std_diff, &max_diff);
		fprintf(data, "%d %d %g %g %g\n", eigenvalue_precision, n_layer_true, mean_diff, std_diff, max_diff);
	}

	eigenvalue_precision = -1;
	for (n_layer = 10; n_layer < n_layer_true; n_layer *= 10)
	{
		n_discontinuity_layer = int(0.25 * n_layer);
		if (4 * n_discontinuity_layer < n_layer) n_discontinuity_layer++;
		printf("N Layer = %i\n", 4 * n_discontinuity_layer);
		probs_test = Precision(eigenvalue_precision, 4 * n_discontinuity_layer);
		Diff(probs_true, probs_test, 1, 1, 1, 0, &mean_diff, &std_diff, &max_diff);
		fprintf(data, "%d %d %g %g %g\n", eigenvalue_precision, 4 * n_discontinuity_layer, mean_diff, std_diff, max_diff);

		n_discontinuity_layer = int(0.25 * 3 * n_layer);
		if (4 * n_discontinuity_layer < 3 * n_layer) n_discontinuity_layer++;
		printf("N Layer = %i\n", 4 * n_discontinuity_layer);
		probs_test = Precision(eigenvalue_precision, 4 * n_discontinuity_layer);
		Diff(probs_true, probs_test, 1, 1, 1, 0, &mean_diff, &std_diff, &max_diff);
		fprintf(data, "%d %d %g %g %g\n", eigenvalue_precision, 4 * n_discontinuity_layer, mean_diff, std_diff, max_diff);
	}

	printf("PREM Full\n");
	n_layer = -1;
	probs_test = Precision(eigenvalue_precision, n_layer);
	Diff(probs_true, probs_test, 1, 1, 1, 0, &mean_diff, &std_diff, &max_diff);
	fprintf(data, "%d %d %g %g %g\n", eigenvalue_precision, n_layer, mean_diff, std_diff, max_diff);

	printf("PREM Four\n");
	n_layer = -2;
	probs_test = Precision(eigenvalue_precision, n_layer);
	Diff(probs_true, probs_test, 1, 1, 1, 0, &mean_diff, &std_diff, &max_diff);
	fprintf(data, "%d %d %g %g %g\n", eigenvalue_precision, n_layer, mean_diff, std_diff, max_diff);

	fclose(data);
}
// We test the impact of the detector depth
void Detector_Depth(int alpha, int beta, bool normal_ordering, bool neutrino_mode, double detector_depth)
{
	double mean_diff, std_diff, max_diff;
	int mo_sign;
	std::vector<std::vector<Matrix3r>> probs_surface, probs_depth;

	// Initialize probability engine
	Probability_Engine probability_engine;

	// Set initial oscillation parameters
	mo_sign = normal_ordering ? +1 : -1;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, mo_sign * 2.534e-3, neutrino_mode); // nu-fit 6

	// Create the energy and cosz vectors
	std::vector<double> Es, coszs;
	double Emin, Emax, Estep, coszmin, coszmax, coszstep;
	int n;

	n = 2e3;

	Emin = 2;
	Emax = 40;

	coszmin = -1.0;
	coszmax = +0.1;

	Estep = (Emax - Emin) / n;
	coszstep = (coszmax - coszmin) / n;

	Es.reserve(n + 1);
	coszs.reserve(n + 1);
	for (int i = 0; i <= n; i++)
	{
		Es.emplace_back(Emin + i * Estep);
		coszs.emplace_back(coszmin + i * coszstep);
	} // i, n

	probability_engine.Set_Spectra(Es, coszs);

	// Create Earth model instance
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);

	// Set Earth details
	// Detector depth
	probability_engine.Set_Production_Height(10);
	probability_engine.Set_Earth(0, &earth_density);
	probs_surface = probability_engine.Get_Probabilities();
	detector_depth = 2;
	probability_engine.Set_Earth(detector_depth, &earth_density);
	probs_depth = probability_engine.Get_Probabilities();

	Diff(probs_surface, probs_depth, alpha, beta, alpha, beta, &mean_diff, &std_diff, &max_diff);
	printf("mean = %g, std = %g, max = %g\n", mean_diff, std_diff, max_diff);

	std::string fname = "data/Detector_Depth_";
	fname += normal_ordering ? "NO" : "IO";
	fname += "_";
	fname += neutrino_mode ? "nu" : "nubar";
	fname += "_";
	fname += (beta == 0) ? "app" : "dis";
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");
	fprintf(data, "%g\n", detector_depth);
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", Es[i]);
	fprintf(data, "\n");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", coszs[i]);
	fprintf(data, "\n");

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
			fprintf(data, "%g ", probs_surface[i][j].arr[alpha][beta] - probs_depth[i][j].arr[alpha][beta]);
		fprintf(data, "\n");
	} // i, n, E

	fclose(data);
}
void Detector_Depth()
{
	double detector_depth = 2;
	int alpha = 1;
	for (int no = 0; no < 2; no++)
	{
		for (int nm = 0; nm < 2; nm++)
		{
			for (int beta = 0; beta < 2; beta++)
				Detector_Depth(alpha, beta, no == 0, nm == 0, detector_depth);
		}
	}
}
// We test the impact of the atmosphere density
void Atmosphere_Density(int alpha, int beta, bool normal_ordering, bool neutrino_mode, double rhoYe_atm)
{
	double mean_diff, std_diff, max_diff;
	int mo_sign;
	std::vector<std::vector<Matrix3r>> probs_vacuum, probs_atmosphere;

	// Initialize probability engine
	Probability_Engine probability_engine;

	// Set initial oscillation parameters
	mo_sign = normal_ordering ? +1 : -1;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, mo_sign * 2.534e-3, neutrino_mode); // nu-fit 6

	// Create the energy and cosz vectors
	std::vector<double> Es, coszs;
	double Emin, Emax, Estep, coszmin, coszmax, coszstep;
	int n;

	n = 2e3;

	Emin = 2;
	Emax = 40;

	coszmin = -1.0;
	coszmax = +0.1;

	Estep = (Emax - Emin) / n;
	coszstep = (coszmax - coszmin) / n;

	Es.reserve(n + 1);
	coszs.reserve(n + 1);
	for (int i = 0; i <= n; i++)
	{
		Es.emplace_back(Emin + i * Estep);
		coszs.emplace_back(coszmin + i * coszstep);
	} // i, n

	probability_engine.Set_Spectra(Es, coszs);

	// Create Earth model instance
	double production_height = 10;
	double detector_depth = 2;
	int n_inner_core_discontinuities = 2;
	int n_outer_core_discontinuities = 10;
	int n_inner_mantle_discontinuities = 10;
	int n_outer_mantle_discontinuities = 5;
	Atmosphere_NDL earth_density_vacuum(n_inner_core_discontinuities, n_outer_core_discontinuities, n_inner_mantle_discontinuities, n_outer_mantle_discontinuities, production_height, 0);
	Atmosphere_NDL earth_density_atmosphere(n_inner_core_discontinuities, n_outer_core_discontinuities, n_inner_mantle_discontinuities, n_outer_mantle_discontinuities, production_height, rhoYe_atm);

	// ** Set Earth details **
	// Vacuum
	probability_engine.Set_Production_Height(0);
	probability_engine.Set_Earth(production_height + detector_depth, &earth_density_vacuum);
	probs_vacuum = probability_engine.Get_Probabilities();

	// Atmosphere
	probability_engine.Set_Earth(production_height + detector_depth, &earth_density_atmosphere);
	probs_atmosphere = probability_engine.Get_Probabilities();

	Diff(probs_vacuum, probs_atmosphere, alpha, beta, alpha, beta, &mean_diff, &std_diff, &max_diff);
	printf("mean = %g, std = %g, max = %g\n", mean_diff, std_diff, max_diff);

	std::string fname = "data/Atmosphere_Density_";
	fname += normal_ordering ? "NO" : "IO";
	fname += "_";
	fname += neutrino_mode ? "nu" : "nubar";
	fname += "_";
	fname += (beta == 0) ? "app" : "dis";
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");
	fprintf(data, "%g\n", rhoYe_atm);
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", Es[i]);
	fprintf(data, "\n");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", coszs[i]);
	fprintf(data, "\n");

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
			fprintf(data, "%g ", probs_vacuum[i][j].arr[alpha][beta] - probs_atmosphere[i][j].arr[alpha][beta]);
		fprintf(data, "\n");
	} // i, n, E

	fclose(data);
}
void Atmosphere_Density()
{
	double rhoYe_atm = 1e-3 * 0.5;
	int alpha = 1;
	for (int no = 0; no < 2; no++)
	{
		for (int nm = 0; nm < 2; nm++)
		{
			for (int beta = 0; beta < 2; beta++)
				Atmosphere_Density(alpha, beta, no == 0, nm == 0, rhoYe_atm);
		}
	}
}
