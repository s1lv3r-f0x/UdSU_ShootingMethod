#define _CRT_SECURE_NO_WARNINGS
#include<cmath>
#include<iostream>

using namespace std;

const int N = 10000;
const double pi = 4 * atan(1);
const double eps = 1e-6, delta = 1e-4;

void func(double x, const double* y, double* dy) {
	const double T = y[0];
	const double q = y[1];
	dy[0] = q;
	dy[1] = -sin(pi * x);
}

int main() {
	double a = 0.0, b = 1.0, dx = 0.0, Tb = 3.0, qa = 1.0;
	double fT = 0.0, t = 0.0, temq = 0.0;
	double T[N + 1], q[N + 1], f0[2], k1[2], k2[2];

	dx = (b - a) / N;
	T[0] = 1.0;
	while (1)
	{
		q[0] = qa;
		for (int i = 0; i < N; i++)
		{
			f0[0] = T[i];
			f0[1] = q[i];
			func(i * dx, f0, k1);
			fT = T[i] + dx * k1[0];
			temq = q[i] + dx * k1[1];
			f0[0] = fT;
			f0[1] = temq;
			func((i + 1) * dx, f0, k2);
			T[i + 1] = T[i] + 0.5 * dx * (k1[0] + k2[0]);
			q[i + 1] = q[i] + 0.5 * dx * (k1[1] + k2[1]);
		}
		//  оличество знаков после зап€той
		cout.precision(10);
		// fabs - модуль числа
		cout << fabs(T[N] - Tb) << "\t|" << T[N] << endl;
		t = T[N];

		if (fabs(T[N] - Tb) < eps) {
			break;
		}

		q[0] += delta;

		for (int i = 0; i < N; i++) {
			f0[0] = T[i];
			f0[1] = q[i];
			func(i * dx, f0, k1);
			fT = T[i] + dx * k1[0];
			temq = q[i] + dx * k1[1];
			f0[0] = fT;
			f0[1] = temq;
			func((i + 1) * dx, f0, k2);
			T[i + 1] = T[i] + 0.5 * dx * (k1[0] + k2[0]);
			q[i + 1] = q[i] + 0.5 * dx * (k1[1] + k2[1]);
		}

		qa -= delta * (t - Tb) / (T[N] - t);
	}
	cout << T[N] << endl;
	system("pause");
	return 0;
}