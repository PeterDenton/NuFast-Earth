#ifndef Earth_H
#define Earth_H

#include <vector>

class Earth_Density
{
	public:
		virtual ~Earth_Density() {};
		int n_discontinuities; // number of discontinuities
		std::vector<double> discontinuities; // the actual discontinuities
		double r_E = 6371; // Earth's radius in km
		bool constant_shells;
		virtual double rhoYe(double r) = 0; // density (g/cc) of the Earth at radius r in km
}; // class Earth_Density

class PREM_NDiscontinuityLayer : public Earth_Density
{
	public:
		PREM_NDiscontinuityLayer(int n_inner_core_discontinuities, int n_outer_core_discontinuities, int n_inner_mantle_discontinuities, int n_outer_mantle_discontinuities);
		PREM_NDiscontinuityLayer(int n_discontinuities);
		double rhoYe(double r);
	private:
		int n_inner_core_discontinuities, n_outer_core_discontinuities, n_inner_mantle_discontinuities, n_outer_mantle_discontinuities;
		double layers[4];
		std::vector<double> rhoYes;
};
class PREM_NUniformLayer : public Earth_Density
{
	public:
		PREM_NUniformLayer(int n_discontinuities);
		double rhoYe(double r);
	private:
		std::vector<double> rhoYes;
};
class PREM_Full : public Earth_Density
{
	public:
		PREM_Full();
		double rhoYe(double r);
};
class PREM_Four : public Earth_Density
{
	public:
		PREM_Four();
		double rhoYe(double r);
};
class PREM_Prob3 : public Earth_Density
{
	public:
		PREM_Prob3();
		double rhoYe(double r);
	private:
		double Ye;
};
class Constant : public Earth_Density
{
	public:
		Constant(double rhoYe_);
		double rhoYe(double r);
	private:
		double rhoYe_;
};

#endif
