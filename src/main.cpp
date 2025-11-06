#include <stdio.h>
#include <vector>
#include <cmath>

#include "NuFastEarth.h"
#include "Matrix.h"

#include "examples/Speed.h"
#include "examples/Figures.h"
#include "examples/Validation.h"

int main()
{
	// Atmospheric neutrino minimal example
	Probability_Engine probability_engine;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, 2.534e-3, true);
	std::vector<double> energies = {1, 2, 3, 4, 5}; // GeV
	std::vector<double> coszs = {-1, -0.5, 0, 1}; // core-crossing to horizontal to down-going
	probability_engine.Set_Spectra(energies, coszs);
	PREM_NDiscontinuityLayer earth_density(2, 10, 10, 5);
	probability_engine.Set_Earth(2, &earth_density); // detector depth in km
	std::vector<std::vector<Matrix3r>> probabilities = probability_engine.Get_Probabilities();
	printf("# E [GeV] | cosz |     Pee     |     Pme     |     Pmm\n");
	for (unsigned int i = 0; i < energies.size(); i++)
	{
		for (unsigned int j = 0; j < coszs.size(); j++)
		{
			printf("%8.2g %6.2g %12.6g %12.6g %12.6g\n", energies[i], coszs[j], probabilities[i][j].arr[0][0], probabilities[i][j].arr[1][0], probabilities[i][j].arr[1][1]);
		} // j, coszs
	} // i, energies

	// Figures
//	Density_Profiles();
//	Oscillogram();
//	Solar_Oscillogram();

	// Validation
//	Eigenvalues_Validation();
//	LBL_Appearance_Validation();
//	Solar_Day_Validation();
//	Solar_Day_Night_Validation();
//	Solar_Weight_Validation();
//	Precision();

	// Speed
//	Geometry_Test_Earth_Density_Model();
//	Eigen_Speed();
//	Atmospheric_Speed();

	return 0;
}
