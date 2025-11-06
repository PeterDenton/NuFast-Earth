#ifndef Geometry_H
#define Geometry_H

#include <vector>

class Earth_Density;

double Distance_Between_Layers(double r_1, double r_2, double sinzsq, double detector_depth, Earth_Density *earth_density);
void Mean_Densities(double cosz, double detector_depth, Earth_Density *earth_density, std::vector<std::pair<double,double>> &mean_densities1, std::vector<std::pair<double,double>> &mean_densities2);

double L2cosz(double L);

#endif
