#pragma once
class AnalitcalSolution {
public:
	// movement
	double u(double x, double t, double l, double So, double ro, double c, int steps);
	// stress 
	double s(double x, double t, double l, double So, double ro, double c,int steps);
	
private:
	double heaviside(double x);
};