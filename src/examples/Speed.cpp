#include <stdio.h>
#include <chrono>
#include <vector>
#include <cmath>
#include <cassert>

#include "examples/Speed.h"
#include "Earth.h"
#include "Geometry.h"
#include "NuFastEarth.h"

#define sq(x) ((x)*(x))

std::string optimization_level = "3"; // change this here as the Makefile is changed

void Geometry_Test_Earth_Density_Model(Earth_Density *ed, int n_1, int n_2)
{
	std::vector<std::pair<double,double>> mean_densities1, mean_densities2; // initialize the distances and densities

	double speed, speed_sum, speedsq_sum;

	speed_sum = 0;
	speedsq_sum = 0;

	for (int i = 0; i < n_1; i++)
	{
		std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < n_2; j++)
			Mean_Densities(-1, 2, ed, mean_densities1, mean_densities2); // calculate the distance and mean density along each section
		std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
		speed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / n_2 * 1e9;
		speed_sum += speed;
		speedsq_sum += sq(speed);
	} // i, n_1
	printf("time = %g +- %g ns\n", speed_sum / n_1, sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1)));
}
void Geometry_Test_Earth_Density_Model()
{
	printf("PREM_Full...\n");
	PREM_Full ed_full;
	Geometry_Test_Earth_Density_Model(&ed_full, 1e3, 1e3);

	printf("PREM_Four...\n");
	PREM_Four ed_four;
	Geometry_Test_Earth_Density_Model(&ed_four, 1e3, 1e3);

	printf("PREM_NUniformLayer(4)...\n");
	PREM_NUniformLayer ed_nuniformlayer0(4);
	Geometry_Test_Earth_Density_Model(&ed_nuniformlayer0, 1e3, 1e3);

	printf("PREM_NUniformLayer(10)...\n");
	PREM_NUniformLayer ed_nuniformlayer1(1e1);
	Geometry_Test_Earth_Density_Model(&ed_nuniformlayer1, 1e3, 1e3);

	printf("PREM_NUniformLayer(100)...\n");
	PREM_NUniformLayer ed_nuniformlayer2(1e2);
	Geometry_Test_Earth_Density_Model(&ed_nuniformlayer2, 1e3, 1e3);

	printf("PREM_NDiscontinuityLayer(2,10,10,5)...\n");
	PREM_NDiscontinuityLayer ed_ndiscontinuitylayer1(2, 10, 10, 5);
	Geometry_Test_Earth_Density_Model(&ed_ndiscontinuitylayer1, 1e3, 1e3);

	printf("PREM_NDiscontinuityLayer(25)...\n");
	PREM_NDiscontinuityLayer ed_ndiscontinuitylayer2(25);
	Geometry_Test_Earth_Density_Model(&ed_ndiscontinuitylayer2, 1e3, 1e3);

	printf("PREM_prob3...\n");
	PREM_Prob3 ed_prob3;
	Geometry_Test_Earth_Density_Model(&ed_prob3, 1e3, 1e3);
}
void Eigen_Speed(int eigenvalue_precision, int n_1, int n_2, FILE *data)
{
	Eigen eigen;
	Probability_Engine probability_engine;

	probability_engine.Set_Eigenvalue_Precision(eigenvalue_precision);
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6
	probability_engine.Precalc();

	double speed, speed_sum, speedsq_sum;

	speed_sum = 0;
	speedsq_sum = 0;

	for (int i = 0; i < n_1; i++)
	{
		std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < n_2; j++)
			eigen = probability_engine.Calculate_Eigen(1.5);
		std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
		speed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / n_2 * 1e9;
		speed_sum += speed;
		speedsq_sum += sq(speed);
	} // i, n_1
	printf("eigenvalue precision = %d, time = %g +- %g ns\n", eigenvalue_precision, speed_sum / n_1, sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1)));
	fprintf(data, "%d %g %g\n", eigenvalue_precision, speed_sum / n_1, sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1)));
} 
void Eigen_Speed()
{
	std::string fname = "data/speed/Eigens_";
	fname += optimization_level;
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");
	for (int eigenvalue_precision = -1; eigenvalue_precision < 3; eigenvalue_precision++)
		Eigen_Speed(eigenvalue_precision, 1e4, 1e4, data);
	fclose(data);
}
// if quick_osc_param: then we cycle through delta or theta23. If not, we cycle through one of the other variables
void Atmospheric_Speed(bool quick_osc_param, int n_1, int n_2, int n_layer, FILE *data)
{
	Probability_Engine probability_engine;
	probability_engine.Set_Oscillation_Parameters(0.307, 0.02195, 0.561, 177 * M_PI / 180, 7.49e-5, +2.534e-3, true); // nu-fit 6

	// Create Earth model instance
	PREM_NUniformLayer earth_density(n_layer);
	// Set Earth details
	probability_engine.Set_Earth(2, &earth_density);
	probability_engine.Set_Production_Height(10);

	// Create the energy and cosz vectors
	std::vector<double> Es = {3};
	std::vector<double> coszs = {-1};
	probability_engine.Set_Spectra(Es, coszs);

	std::vector<std::vector<Matrix3r>> probs;
	double speed, speed_sum, speedsq_sum;

	speed_sum = 0;
	speedsq_sum = 0;

	for (int i = 0; i < n_1; i++)
	{
		std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < n_2; j++)
		{
			if (quick_osc_param)
				probability_engine.Set_s23sq(0.45 + i * 0.1 / n_2);
			else
				probability_engine.Set_Dmsq31(2.5e-3 + i * 0.1e-3 / n_2);
			probs = probability_engine.Get_Probabilities();
		} // j, n_2
		std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
		speed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / n_2 * 1e9;
		speed_sum += speed;
		speedsq_sum += sq(speed);
	} // i, n_1
	if (quick_osc_param)
	{
		fprintf(data, "quick ");
		printf("quick ");
	}
	else
	{
		fprintf(data, "full ");
		printf("full ");
	}
	printf("n_layer = %d, time = %g +- %g ns\n", n_layer, speed_sum / n_1, sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1)));
	fprintf(data, "%d %g %g\n", n_layer, speed_sum / n_1, sqrt(speedsq_sum / n_1 - sq(speed_sum / n_1)));
}
void Atmospheric_Speed()
{
	const int n_n_layer = 4;
	int n_layers[n_n_layer] = {4, 10, 100, 1000};

	std::string fname = "data/speed/Atmospherics_";
	fname += optimization_level;
	fname += ".txt";
	FILE *data = fopen(fname.c_str(), "w");
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < n_n_layer; j++)
			Atmospheric_Speed(i == 0, 1e4, 1e2, n_layers[j], data);
	} // quick/slow
	fclose(data);
}
