#include <cassert>
#include <cmath>
#include <vector>

#include "Geometry.h"
#include "Earth.h"

#define sq(x) ((x)*(x))

double Distance_Between_Layers(double r_1, double r_2, double sinzsq, double detector_depth, Earth_Density *earth_density)
{
	assert(r_1 > r_2);

	double a, bsq, csq;

    a = earth_density->r_E - detector_depth;
    bsq = sq(r_1 / a) - sinzsq;
    csq = sq(r_2 / a) - sinzsq;
    if (bsq < 0 and bsq > -1e-14) bsq = 0;
    if (csq < 0 and csq > -1e-14) csq = 0;
    return a * (sqrt(bsq) - sqrt(csq));
}

// given r1, what radius is L farther away?
double L2r(double L, double r_1, double sinzsq, double detector_depth, Earth_Density *earth_density)
{
    return sqrt(sq(L) + sq(r_1) - 2 * L * (earth_density->r_E - detector_depth) * sqrt(sq(r_1) / sq(earth_density->r_E - detector_depth) - sinzsq));
}

void Mean_Density(double sinzsq, double detector_depth, double r_1, double r_f, Earth_Density *earth_density, std::vector<std::pair<double,double>> &mean_density)
{
	double r, r_2, L, L_step, L_step_, rhoYe_sum;
	int n_step, i_discontinuity;

	L_step = 10.; // approximate step size in km for the integral

	mean_density.clear();
	mean_density.reserve(2 * earth_density->n_discontinuities + 2);

	while (r_1 > r_f)
	{
		r_2 = r_f;
		// check if there is a discontinuity to hit before the target
		i_discontinuity = -1;
		for (int i = earth_density->n_discontinuities - 1; i >= 0; i--)
		{
			if (earth_density->discontinuities[i] < r_1)
			{
				if (i_discontinuity == -1) i_discontinuity = i;
				if (earth_density->discontinuities[i] > r_f)
				{
					r_2 = earth_density->discontinuities[i];
					break;
				}
			}
		} // i_discontinuity, n_discontinuities, discontinuity
		L = Distance_Between_Layers(r_1, r_2, sinzsq, detector_depth, earth_density);
		// set up integral
		if (earth_density->constant_shells)
		{
			mean_density.emplace_back(L, i_discontinuity + 1);
		} // density does not vary in a shell, do not do an integral
		else
		{
			n_step = int(ceil(1. + L / L_step));
			L_step_ = L / n_step; // actual step size with rounding handled
			// do integral
			rhoYe_sum = 0;
			for (int i = 0; i < n_step; i++)
			{
				r = L2r((i + 0.5) * L_step_, r_1, sinzsq, detector_depth, earth_density);
				rhoYe_sum += earth_density->rhoYe(r);
			} // i, n_step, L
			mean_density.emplace_back(L, rhoYe_sum / n_step);
		} // density does vary within a shel, do an integral
		r_1 = r_2;
	}
}

// cosz: cos(theta_z) which is -1 for core crossing and +1 for coming from above
// detector depth (km) is zero for a detector at the surface and is positive if the detector is in the Earth
void Mean_Densities(double cosz, double detector_depth, Earth_Density *earth_density, std::vector<std::pair<double,double>> &mean_densities1, std::vector<std::pair<double,double>> &mean_densities2)
{
	assert(sq(cosz) <= 1);
	assert(detector_depth >= 0);

	double sinzsq, sinz, r_min, r_1, r_f;

	sinzsq = 1 - sq(cosz);
	sinz = sqrt(sinzsq);
	r_min = (earth_density->r_E - detector_depth) * sinz;
	if (cosz > 0) r_min = earth_density->r_E - detector_depth; // handles down going trajectories and only includes the first set of mean_densities

	// go from the surface to the detector depth
	r_1 = earth_density->r_E; // always assume we begin at the earth's surface: add atmosphere later in the calculation
	r_f = earth_density->r_E - detector_depth; // end

	Mean_Density(sinzsq, detector_depth, r_1, r_f, earth_density, mean_densities1);

	// go from detector depth to the halfway point
	r_1 = r_f;
	r_f = r_min;
	Mean_Density(sinzsq, detector_depth, r_1, r_f, earth_density, mean_densities2);
}

// converts a distance for a LBL experiment (km) to cosz
double L2cosz(double L)
{
	assert (L >= 0);
	assert (L <= 2 * 6371.);

	return -L / (2 * 6371.);
}
