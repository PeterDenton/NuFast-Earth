#include <cmath>
#include <cassert>

#include "Solar.h"

#define sq(x) ((x)*(x))

// AGM formula for the complete elliptic integral to eliminate
// the dependency on std::comp_ellint_1 from cmath which is not
// available on macOS (yet)
static double comp_ellint_1(double k)
{
	assert (k <= 1);

	double a = 1.0;
	double b = std::sqrt(1.0 - sq(k));
	double c;

	while (fabs(a - b) > 1e-15) // 1e-16 is too small and it may not converge
	{
		c = (a + b) * 0.5;
		b = std::sqrt(a * b);
		a = c;
	}
	return M_PI / (2.0 * a);
}


// table II in hep-ph/9702343
// eta is the nadir angle in radians
// latitude is the experiments latitude in radians
double Solar_Weight(double eta, double latitude)
{
	if (eta > 0.5 * M_PI) return 0; // zero for nadirs above the horizon
	double i, y, z, num, den;

	i = 0.4091; // Earth's inclination, see after eq. c4
/*  latitude = 46.47 * M_PI / 180; // SNO latitude, see fig. 3
	latitude = 44.356 * M_PI / 180; // DUNE far detector latitude, see end of table 1 in 1707.02322
	latitude = 42.45 * M_PI / 180; // Gran Sasso (Borexino) latitude, see fig. 3 hep-ph/9702343
	latitude = 36.48 * M_PI / 180; // kamioka, see fig. 3 hep-ph/9702343
*/
	if (eta < latitude - i) // first row
		return 0;

	z = sqrt(sin(i) * cos(latitude) * sin(eta));
	y = sqrt(
			sin(0.5 * (i + latitude + eta)) * 
			sin(0.5 * (i - latitude + eta)) * 
			cos(0.5 * (i + latitude - eta)) * 
			cos(0.5 * (i - latitude - eta)));


	if (latitude < i) // first column
	{
		if (i - latitude < eta and eta < i + latitude)
		{
			num = y;
			den = z;
		}
		else
		{
			num = z;
			den = y;
		}
	} // first column
	else if (latitude < 0.5 * M_PI - i) // second column
	{
		if (eta > latitude + i)
		{
			num = z;
			den = y;
		}
		else
		{
			num = y;
			den = z;
		}
	} // second column
	else // third column
	{
		if (eta > M_PI - latitude - i)
		{
			num = z;
			den = y;
		}
		else
		{
			num = y;
			den = z;
		}
	} // third column

	return 2 * sin(eta) * comp_ellint_1(num / den) / (sq(M_PI) * den);
}
