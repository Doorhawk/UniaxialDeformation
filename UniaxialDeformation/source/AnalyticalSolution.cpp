#include "AnalyticalSolution.h"

double AnalitcalSolution::u(double x, double t, double l, double So, double ro, double c, int steps){
	double answ = 0;
	int sign = 1;
	for (int i = 0; i < steps; i++, sign *= -1) {
		double xminus = t - ((2 * i + 1) * l - x) / c;
		double xplus = t - ((2 * i + 1) * l + x) / c;
		answ += sign * (xminus * heaviside(xminus) - xplus * heaviside(xplus));
	}

	return answ * (-So / ro / c);
}
double AnalitcalSolution::s(double x, double t, double l, double So, double ro, double c, int steps){
	double answ = 0;
	int sign = 1;
	for (int i = 0; i < steps; i++, sign *= -1) {
		double xminus = t - ((2 * i + 1) * l - x) / c;
		double xplus = t - ((2 * i + 1) * l + x) / c;
		answ += sign * (heaviside(xminus) + heaviside(xplus));
	}

	return answ * (-So);
}
double AnalitcalSolution::heaviside(double x) {
	return x >= 0 ? 1 : 0;
}
