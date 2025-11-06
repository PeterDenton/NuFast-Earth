#ifndef Matrix_H
#define Matrix_H

#include <complex>

#define sq(x) ((x)*(x))

class Matrix3r
{
	public:
		Matrix3r()
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					arr[i][j] = 0;
			} // i
		};
		~Matrix3r() {};
		double arr[3][3];
		void Print(int n)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					printf("%.16e ", arr[i][j]);
				printf("\n");
			} // i
			for (int i = 0; i < n; i++)
				printf("\n");
		};
		void Print() { Print(0); };

		Matrix3r operator*(Matrix3r const& m)
		{
			Matrix3r ret;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
						ret.arr[i][j] += arr[i][k] * m.arr[k][j];
				} // j
			} // i
			return ret;
		};
		Matrix3r Transpose()
		{
			Matrix3r ret;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					ret.arr[i][j] = arr[j][i];
			} // i
			return ret;
		}
};

class Vector3r
{
	public:
		Vector3r()
		{
			for (int i = 0; i < 3; i++)
				vec[i] = 0;
		};
		~Vector3r() {};
		double vec[3];
		void Print(int n)
		{
			for (int i = 0; i < 3; i++)
				printf("%g ", vec[i]);
			for (int i = 0; i < n; i++)
				printf("\n");
		};
		void Print() { Print(0); };
};

class Matrix3c
{
	public:
		Matrix3c()
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					arr[i][j] = 0;
			} // i
		};
		~Matrix3c() {};
		std::complex<double> arr[3][3];
		void Print(int n)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					printf("%.16f+%.16fi ", arr[i][j].real(), arr[i][j].imag());
				printf("\n");
			} // i
			for (int i = 0; i < n; i++)
				printf("\n");
		};
		void Print() { Print(0); };

		void Identity()
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
					arr[i][j] = 0;
			} // i
			for (int i = 0; i < 3; i++)
				arr[i][i] = 1;
		}
		Matrix3c operator*(Matrix3c const& m)
		{
			Matrix3c ret;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
						ret.arr[i][j] += arr[i][k] * m.arr[k][j];
				} // j
			} // i
			return ret;
		};

		// returns the matrix times the transpose of the matrix
		Matrix3c AAT()
		{
			Matrix3c ret;
			for (int i = 0; i < 3; i++)
			{
				ret.arr[i][i] = sq(arr[i][0]) + sq(arr[i][1]) + sq(arr[i][2]);
				for (int j = i + 1; j < 3; j++)
				{
					ret.arr[i][j] = arr[i][0] * arr[j][0] + arr[i][1] * arr[j][1] + arr[i][2] * arr[j][2];
					ret.arr[j][i] = ret.arr[i][j];
				} // j
			} // i, 3
			return ret;
		};
};

struct Eigen
{
	Vector3r lambda;
	double Ve2sq, Ve3sq, Vm2sq, Vm3sq, Vt2sq, Vt3sq, Ve3Vm3, Ve3Vt3, Vm3Vt3, Ve2Vm2, Ve2Vt2, Vm2Vt2;
};

#endif
