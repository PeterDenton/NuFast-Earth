#include <stdio.h>
#include <vector>

#include "examples/Figures.h"
#include "Earth.h"
#include "Matrix.h"
#include "NuFastEarth.h"

void Density_Profiles()
{
	PREM_Full earth_density_full;
	PREM_NDiscontinuityLayer earth_density_ndiscontinuitylayer(2, 10, 10, 5);
	int n_layer = 20;
	PREM_NUniformLayer earth_density_nuniformlayer(n_layer);
	PREM_Four earth_density_four;
	PREM_Prob3 earth_density_prob3;

	double r, rmin, rmax, rstep;
	int n = 1e4;

	rmin = 0;
	rmax = 6371;
	rstep = (rmax - rmin) / n;

	FILE *data = fopen("data/Density_Profiles.txt", "w");
	fprintf(data, "%d\n", n_layer);
	for (int i = 0; i <= n; i++)
	{
		r = rmin + i * rstep;
		fprintf(data, "%g ", r);
		fprintf(data, "%g ", earth_density_full.rhoYe(r));
		fprintf(data, "%g ", earth_density_ndiscontinuitylayer.rhoYe(r));
		fprintf(data, "%g ", earth_density_nuniformlayer.rhoYe(r));
		fprintf(data, "%g ", earth_density_four.rhoYe(r));
		fprintf(data, "%g ", earth_density_prob3.rhoYe(r));
		fprintf(data, "\n");
	} // i, n, r
	fclose(data);
}
void Oscillogram(int alpha, int beta, bool neutrino_mode, bool normal_ordering)
{
	Probability_Engine probability_engine;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, neutrino_mode); // nu-fit 6

	if (not normal_ordering)	probability_engine.Set_Dmsq31(-2.534e-3);

	// Create Earth model instance
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);
	// Set Earth details
	probability_engine.Set_Earth(2, &earth_density);
	probability_engine.Set_Production_Height(10);

	// Create the energy and cosz vectors
	std::vector<double> Es, coszs;
	double Emin, Emax, Estep, coszmin, coszmax, coszstep;
	int n;

	n = 1e3;

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

	std::vector<std::vector<Matrix3r>> probs;
	probs = probability_engine.Get_Probabilities();

	std::string fname = "data/Oscillogram_";
	fname += std::to_string(alpha) + "_" + std::to_string(beta) + "_";
	fname += neutrino_mode ? "nu" : "nubar";
	fname += "_";
	fname += normal_ordering ? "NO" : "IO";
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", Es[i]);
	fprintf(data, "\n");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", coszs[i]);
	fprintf(data, "\n");

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
			fprintf(data, "%g ", probs[i][j].arr[alpha][beta]);
		fprintf(data, "\n");
	} // i, n, E

	fclose(data);
}
void Oscillogram()
{
	for (int beta = 0; beta < 2; beta++)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				Oscillogram(1, beta, i == 0, j == 0);
		} // i, 2, neutrino mode
	}
}

void Solar_Oscillogram()
{
	Probability_Engine probability_engine;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6

	// Create Earth model instance
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);
	// Set Earth details
	probability_engine.Set_Earth(2, &earth_density);

	// Create the energy and cosz vectors
	std::vector<double> Es, coszs;
	double Emin, Emax, Estep, coszmin, coszmax, coszstep;
	int n;

	n = 1e3;

	Emin = 5e-3;
	Emax = 20e-3;

	coszmin = -1;
	coszmax = +0.2;

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

	probability_engine.Set_rhoYe_Sun(100 * (2. / 3));

	std::vector<std::vector<Matrix3r>> probs;
	probs = probability_engine.Get_Solar_Night_Probabilities();

	FILE *data = fopen("data/Solar_Oscillogram.txt", "w");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", Es[i]);
	fprintf(data, "\n");
	for (int i = 0; i <= n; i++)
		fprintf(data, "%g ", coszs[i]);
	fprintf(data, "\n");

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
			fprintf(data, "%g ", probs[i][j].arr[0][0]);
		fprintf(data, "\n");
	} // i, n, E

	fclose(data);
}
