#include <cassert>
#include <cmath>
#include <vector>

#include "Earth.h"

#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

double PREM_Full_rho(double r)
{
	double x;

	x = r / 6371.;

	if (r < 1221.5) return 13.0885 - 8.8381 * sq(x);
	if (r < 3480.0) return 12.5815 - 1.2638 * x - 3.6426 * sq(x) - 5.5281 * cube(x);
	if (r < 5701.0) return 7.9565 - 6.4761 * x + 5.5283 * sq(x) - 3.0807 * cube(x);
	if (r < 5771.0) return 5.3197 - 1.4836 * x;
	if (r < 5971.0) return 11.2494 - 8.0298 * x;
	if (r < 6151.0) return 7.1089 - 3.8045 * x;
	if (r < 6346.6) return 2.6910 + 0.6924 * x;
	if (r < 6356.0) return 2.9;
	if (r < 6368.0) return 2.6;
	if (r <= 6371.0) return 1.02;
	return 0;
}
double PREM_Full_Ye(double r)
{
	// Set Ye from hep-ph/0002149 above table 1
	if (r < 3480.0)	return 0.466;
	else			return 0.494;
}
///////////////////////////////////////////////////////////////////////////////////////////
PREM_NDiscontinuityLayer::PREM_NDiscontinuityLayer(int n_inner_core_discontinuities_, int n_outer_core_discontinuities_, int n_inner_mantle_discontinuities_, int n_outer_mantle_discontinuities_)
{
	n_inner_core_discontinuities = n_inner_core_discontinuities_;
	n_outer_core_discontinuities = n_outer_core_discontinuities_;
	n_inner_mantle_discontinuities = n_inner_mantle_discontinuities_;
	n_outer_mantle_discontinuities = n_outer_mantle_discontinuities_;

	n_discontinuities = n_inner_core_discontinuities + n_outer_core_discontinuities + n_inner_mantle_discontinuities + n_outer_mantle_discontinuities;
	discontinuities.reserve(n_discontinuities);
	rhoYes.reserve(n_discontinuities);

	layers[0] = 1221.5;
	layers[1] = 3480;
	layers[2] = 5701;
	layers[3] = 6371;

	double r;
	int count;

	count = 0;
	// inner core
	for (int i = 0; i < n_inner_core_discontinuities; i++)
	{
		discontinuities[count] = (i + 1) * layers[0] / n_inner_core_discontinuities;
		r = (i + 0.5) * layers[0] / n_inner_core_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_inner_core_discontinuities, r
	// outer core
	for (int i = 0; i < n_outer_core_discontinuities; i++)
	{
		discontinuities[count] = layers[0] + (i + 1) * (layers[1] - layers[0]) / n_outer_core_discontinuities;
		r = layers[0] + (i + 0.5) * (layers[1] - layers[0]) / n_outer_core_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_outer_core_discontinuities, r
	// inner mantle
	for (int i = 0; i < n_inner_mantle_discontinuities; i++)
	{
		discontinuities[count] = layers[1] + (i + 1) * (layers[2] - layers[1]) / n_inner_mantle_discontinuities;
		r = layers[1] + (i + 0.5) * (layers[2] - layers[1]) / n_inner_mantle_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_inner_mantle_discontinuities, r
	// outer mantle
	for (int i = 0; i < n_outer_mantle_discontinuities; i++)
	{
		discontinuities[count] = layers[2] + (i + 1) * (layers[3] - layers[2]) / n_outer_mantle_discontinuities;
		r = layers[2] + (i + 0.5) * (layers[3] - layers[2]) / n_outer_mantle_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_outer_mantle_discontinuities, r

	constant_shells = true;
}
PREM_NDiscontinuityLayer::PREM_NDiscontinuityLayer(int n_discontinuities_) : PREM_NDiscontinuityLayer(n_discontinuities_, n_discontinuities_, n_discontinuities_, n_discontinuities_) {}
double PREM_NDiscontinuityLayer::rhoYe(double r)
{
	assert(r >= 0);
	if (r > layers[3]) return 0;
	if (r > layers[2])
		return rhoYes[(int)floor((r - layers[2]) / (layers[3] - layers[2]) * n_outer_mantle_discontinuities) + n_inner_core_discontinuities + n_outer_core_discontinuities + n_inner_mantle_discontinuities];
	if (r > layers[1])
		return rhoYes[(int)floor((r - layers[1]) / (layers[2] - layers[1]) * n_inner_mantle_discontinuities) + n_inner_core_discontinuities + n_outer_core_discontinuities];
	if (r > layers[0])
		return rhoYes[(int)floor((r - layers[0]) / (layers[1] - layers[0]) * n_inner_mantle_discontinuities) + n_inner_core_discontinuities];
	return rhoYes[(int)floor(r / layers[0] * n_inner_core_discontinuities)];
}
///////////////////////////////////////////////////////////////////////////////////////////
PREM_NUniformLayer::PREM_NUniformLayer(int n_discontinuities_)
{
	n_discontinuities = n_discontinuities_;
	discontinuities.reserve(n_discontinuities);
	rhoYes.reserve(n_discontinuities);

	double r;

	for (int i = 0; i < n_discontinuities; i++)
	{
		discontinuities[i] = (i + 1) * 6371. / n_discontinuities;
		r = (i + 0.5) * 6371. / n_discontinuities;
		rhoYes[i] = PREM_Full_Ye(r) * PREM_Full_rho(r);
	} // i, n_discontinuities, r

	constant_shells = true;
}
double PREM_NUniformLayer::rhoYe(double r)
{
	assert(r >= 0);
	if (r > 6371.) return 0;
	return rhoYes[(int)floor(r / r_E * n_discontinuities)];
}
///////////////////////////////////////////////////////////////////////////////////////////
PREM_Full::PREM_Full()
{
	n_discontinuities = 10;
	discontinuities.reserve(n_discontinuities);
	discontinuities[0] = 1221.5;
	discontinuities[1] = 3480.;
	discontinuities[2] = 5701.;
	discontinuities[3] = 5771.;
	discontinuities[4] = 5971.;
	discontinuities[5] = 6151.;
	discontinuities[6] = 6346.6;
	discontinuities[7] = 6356.;
	discontinuities[8] = 6368.;
	discontinuities[9] = 6371.;

	constant_shells = false;
}
double PREM_Full::rhoYe(double r)
{
	assert(r >= 0);

	return PREM_Full_Ye(r) * PREM_Full_rho(r);
}
///////////////////////////////////////////////////////////////////////////////////////////
PREM_Four::PREM_Four()
{
	n_discontinuities = 4;
	discontinuities.reserve(n_discontinuities);
	discontinuities[0] = 1221.5;
	discontinuities[1] = 3480.;
	discontinuities[2] = 5701.;
	discontinuities[3] = 6371.;

	constant_shells = false;
}
double PREM_Four::rhoYe(double r)
{
	assert(r >= 0);

	double x;

	x = r / r_E;

	if (r < 1221.5) return 0.466 * (13.0885 - 8.8381 * sq(x));
	if (r < 3480.0) return 0.466 * (12.5815 - 1.2638 * x - 3.6426 * sq(x) - 5.5281 * cube(x));
	if (r < 5701.0) return 0.494 * (7.9565 - 6.4761 * x + 5.5283 * sq(x) - 3.0807 * cube(x));
	if (r <= 6371.0) return 0.494 * (5.3197 - 1.4836 * x);
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////
PREM_Prob3::PREM_Prob3()
{
	n_discontinuities = 4;
	discontinuities.reserve(n_discontinuities);
	discontinuities[0] = 1220.;
	discontinuities[1] = 3480.;
	discontinuities[2] = 5701.;
	discontinuities[3] = 6371.;

	Ye = 0.5;

	constant_shells = true;
}
double PREM_Prob3::rhoYe(double r)
{
	assert(r >= 0);

	if (r > 6371.) return 0.;
	if (r > 5701.) return Ye * 3.3;
	if (r > 3480.) return Ye * 5.;
	if (r > 1220.) return Ye * 11.3;
	return Ye * 13.;
}
///////////////////////////////////////////////////////////////////////////////////////////
Constant::Constant(double rhoYe__)
{
	n_discontinuities = 1;
	discontinuities.reserve(n_discontinuities);
	discontinuities[0] = 6371.;

	constant_shells = true;

	rhoYe_ = rhoYe__;
}
double Constant::rhoYe(double r)
{
	assert(r >= 0);

    if (r <= 6371.) return rhoYe_;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////
// production_height in NuFastEarth object should be zero
Atmosphere_NDL::Atmosphere_NDL(int n_inner_core_discontinuities_, int n_outer_core_discontinuities_, int n_inner_mantle_discontinuities_, int n_outer_mantle_discontinuities_, double production_height_, double rhoYe_atm_)
{
	n_inner_core_discontinuities = n_inner_core_discontinuities_;
	n_outer_core_discontinuities = n_outer_core_discontinuities_;
	n_inner_mantle_discontinuities = n_inner_mantle_discontinuities_;
	n_outer_mantle_discontinuities = n_outer_mantle_discontinuities_;
	production_height = production_height_;
	rhoYe_atm = rhoYe_atm_;

	r_E = 6371 + production_height;

	n_discontinuities = n_inner_core_discontinuities + n_outer_core_discontinuities + n_inner_mantle_discontinuities + n_outer_mantle_discontinuities + 1;
	discontinuities.reserve(n_discontinuities);
	rhoYes.reserve(n_discontinuities);

	layers[0] = 1221.5;
	layers[1] = 3480;
	layers[2] = 5701;
	layers[3] = 6371;
	layers[4] = 6371 + production_height;

	double r;
	int count;

	count = 0;
	// inner core
	for (int i = 0; i < n_inner_core_discontinuities; i++)
	{
		discontinuities[count] = (i + 1) * layers[0] / n_inner_core_discontinuities;
		r = (i + 0.5) * layers[0] / n_inner_core_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_inner_core_discontinuities, r
	// outer core
	for (int i = 0; i < n_outer_core_discontinuities; i++)
	{
		discontinuities[count] = layers[0] + (i + 1) * (layers[1] - layers[0]) / n_outer_core_discontinuities;
		r = layers[0] + (i + 0.5) * (layers[1] - layers[0]) / n_outer_core_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_outer_core_discontinuities, r
	// inner mantle
	for (int i = 0; i < n_inner_mantle_discontinuities; i++)
	{
		discontinuities[count] = layers[1] + (i + 1) * (layers[2] - layers[1]) / n_inner_mantle_discontinuities;
		r = layers[1] + (i + 0.5) * (layers[2] - layers[1]) / n_inner_mantle_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_inner_mantle_discontinuities, r
	// outer mantle
	for (int i = 0; i < n_outer_mantle_discontinuities; i++)
	{
		discontinuities[count] = layers[2] + (i + 1) * (layers[3] - layers[2]) / n_outer_mantle_discontinuities;
		r = layers[2] + (i + 0.5) * (layers[3] - layers[2]) / n_outer_mantle_discontinuities;
		rhoYes[count] = PREM_Full_Ye(r) * PREM_Full_rho(r);
		count++;
	} // i, n_outer_mantle_discontinuities, r
	
	discontinuities[count] = layers[4];
	rhoYes[count] = rhoYe_atm; // atmosphere density is about 1e-3 g/cm^3, and the electron fraction of both nitrogen and oxygen is 0.5

	constant_shells = true;
}
Atmosphere_NDL::Atmosphere_NDL(int n_discontinuities_, double production_height_, double rhoYe_atm_) : Atmosphere_NDL(n_discontinuities_, n_discontinuities_, n_discontinuities_, n_discontinuities_, production_height_, rhoYe_atm_) {}
double Atmosphere_NDL::rhoYe(double r)
{
	assert(r >= 0);
	if (r > layers[4]) assert(false);
	if (r > layers[3])
		return rhoYes[n_inner_core_discontinuities + n_outer_core_discontinuities + n_inner_mantle_discontinuities + n_outer_mantle_discontinuities];
	if (r > layers[2])
		return rhoYes[(int)floor((r - layers[2]) / (layers[3] - layers[2]) * n_outer_mantle_discontinuities) + n_inner_core_discontinuities + n_outer_core_discontinuities + n_inner_mantle_discontinuities];
	if (r > layers[1])
		return rhoYes[(int)floor((r - layers[1]) / (layers[2] - layers[1]) * n_inner_mantle_discontinuities) + n_inner_core_discontinuities + n_outer_core_discontinuities];
	if (r > layers[0])
		return rhoYes[(int)floor((r - layers[0]) / (layers[1] - layers[0]) * n_inner_mantle_discontinuities) + n_inner_core_discontinuities];
	return rhoYes[(int)floor(r / layers[0] * n_inner_core_discontinuities)];
}

