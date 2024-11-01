#pragma once
#include "AnalyticalSolution.h"
#include "comma.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <iostream>

using namespace std;

class FEM {
public:
    FEM(int _n, double _l, double _ro0, double _K, double _G, double _Tmax);

    void set_border_v(bool left1, double v0_l1, bool right1, double v0_r1);
    void set_border_f(bool left1, double f0_l1, bool right1, double f0_r1);
    void set_print_parametres(bool printAnalitcalSolution, bool numberFormatIsdot, string& separator);

    void solve();

private:
    int n;                          // n - number of nodes, (n-1) - number of cells
    double Tmax;                    // Tmax - end time
    double K;                       // K - bulk modulus
    double G;                       // G - shear modulus
    double E;                       // E - Young's modulus
    double Nu;                      // Nu - Poisson's ratio
    double c;                       // c - speed of sound
    double l;                       // l - length
    double ro0;                     // ro0 - initial density
    double t;                       // t - time
    double dt;                      // dt - time step distance between ti+1 - ti
    double dt2;                     // dt2 - half time step distance between ti+1/2 - ti-1/2
    vector<double> m_cell;          // m_cell - cell mass;
    vector<double> m_dot;           // m_dot - nodal mass
    vector<double> ro;              // ro - actual cell density
    vector<double> x0;              // x0 - initial node coordinates
    vector<double> x;               // x - node coordinates with time step dt
    vector<double> x2;              // x2 - node coordinates with time step dt2
    vector<double> v;               // v - node velocity
    vector<double> u;               // u - node displacement
    vector<double> f;               // f - node forces 
    vector<double> psigm;           // psigm - spherical stress tensor
    vector<vector<double>> dsigm;   // dsigm - stress tensor deviator
    vector<vector<double>> sigma;   // sigma - full stress tensor


    bool is_v0_left, is_v0_right;   // are there boundary conditions on the left/right end for VELOCITY
    double v0_l, v0_r;              // if there are, what is the speed set on it

    bool is_f0_left, is_f0_right;   // are there boundary conditions on the left/right end for FORCE
    double f0_l, f0_r;              // if there are, what is the force set on it

    bool printAS;                   // print analytical solution
    bool numberFormat;              // print number in format 0.00 - true 0,00 - false
    string sep;                     // separator in output

    AnalitcalSolution ana;
    fstream ffu;
    fstream ffs;
    int printCounter;

    void printSelectedInformation();
    void step();

    void set_border_v();
    void set_border_f();
    void set_start();

    void calc_dti_ro();
    void calc_x_u();
    void calc_eps_sigm();
    void calc_forse();
    void calc_v();
    void correct_v();
};