#include "FiniteElementMethod.h"



void FEM::set_border_v() {
    if (is_v0_left) {
        v[0] = v0_l;
    }
    if (is_v0_right) {
        v[n - 1] = v0_r;
    }
}
void FEM::set_border_f() {
    if (is_f0_left) {
        f[0] = f[0] + f0_l;
    }
    if (is_f0_right) {
        f[n - 1] = f[n - 1] + f0_r;
    }
}
void FEM::set_start() {
}

void FEM::calc_dti_ro() {

    for (int i = 0; i < n-1; i++) {
        ro[i] = m_cell[i] / (x[i + 1] - x[i]) * 2;
    }

    double min_t = FLT_MAX, t1 = 0;
    for (int i = 0; i < n - 1; i++) {
        t1 = sqrt(ro[i] / (K + 4. / 3 * G)) * (x[i + 1] - x[i]);
        if (t1 < min_t) {
            min_t = t1;
        }
    }

    double tmp = dt;
    dt = min_t * 0.9;
    dt2 = (dt + tmp) / 2;

    t += dt;
}
void FEM::calc_x_u() {
    double u2;
    for (int i = 0; i < n; i++) {
        u2 = u[i] + dt * v[i] / 2;
        u[i] = u[i] + dt * v[i];
        x2[i] = x0[i] + u2;
        x[i] = x0[i] + u[i];
    }
}
void FEM::calc_eps_sigm() {
    double eps0v;
    for (int i = 0; i < n - 1; i++) {
        eps0v = (v[i + 1] - v[i]) / (x2[i + 1] - x2[i]);

        psigm[i] = psigm[i] - K * eps0v * dt2;     

        dsigm[i][0] = dsigm[i][0] + 2 * G * (2. / 3 * eps0v) * dt2;                 
        dsigm[i][1] = dsigm[i][2] = dsigm[i][1] - 2 * G * (1. / 3 * eps0v) * dt2;

        for (int j = 0; j < 3; j++) {
            sigma[i][j] = dsigm[i][j] - psigm[i];                                   
        }
    }
}
void FEM::calc_forse() {
    f[0] = -sigma[0][0];
    for (int i = 1; i < n - 1; i++)
        f[i] = sigma[i - 1][0] - sigma[i][0];
    f[n - 1] = sigma[n - 2][0];
}
void FEM::calc_v() {
    for (int i = 0; i < n; i++)
        v[i] = (v[i] - f[i] * dt2 / m_dot[i]);
}
void FEM::correct_v() {
    vector<double> fi2(n - 1);
    vector<double> fi(n);

    for (int i = 0; i < n - 1; i++) {
        fi2[i] = 1. / 12 * m_cell[i] * (v[i + 1] - v[i]);
    }

    fi[0] = fi2[0];
    for (int i = 1; i < n - 1; i++) {
        fi[i] = fi2[i] - fi2[i - 1];
    }
    fi[n - 1] = fi2[n - 2];

    for (int i = 0; i < n; i++) {
        v[i] = v[i] + fi[i] / m_dot[0];
    }
}


void FEM::step() {

    calc_dti_ro();
    calc_x_u();
    calc_eps_sigm();
    calc_forse();
    set_border_f();
    calc_v();
    set_border_v();
    correct_v();
    set_border_v();

}
void FEM::printSelectedInformation() {

    cout << t << endl;

    ffs << t << ";";
    ffs << sigma[0][0] << ";";
    ffs << sigma[(n - 2) / 2][0] << ";";
    if (printAS) {
        ffs << ana.s(0, t, l, f0_r, ro0, c,10000) << ";";
        ffs << ana.s(l / 2, t, l, f0_r, ro0, c,10000) << ";";
    }
    ffs << endl;


    ffu << t << ";";
    ffu << u[n - 1] << ";";
    ffu << u[(n - 1) / 2] << ";";
    if (printAS) {
        ffu << ana.u(l, t, l, f0_r, ro0, c,10000) << ";";
        ffu << ana.u(l / 2, t, l, f0_r, ro0, c,10000) << ";";
    }
    ffu << endl;

}


void FEM::solve(double Tmax) {
    while (t < Tmax) {
        printSelectedInformation();
        step();
    }
}

FEM::FEM(int _n, double _l, double _ro0, double _K, double _G) {

    n = _n;
    l = _l;
    ro0 = _ro0;
    G = _G;   
    K = _K;

    E = 9 * K * G / (3 * K + G);
    Nu = (3 * K - 2 * G) / (2 * (3 * K + G));
    c = sqrt(E / ro0 * (1 - Nu) / (1 + Nu) / (1 - 2 * Nu));


    x0.resize(n);
    x.resize(n);
    x2.resize(n);
    v.resize(n);
    u.resize(n);
    f.resize(n);
    for (int i = 0; i < n; i++) {
        x0[i] = x[i] = double(i) / (n - 1) * l;  
    }

    psigm.resize(n-1);
    m_cell.resize(n - 1);
    ro.resize(n - 1);
    dsigm.resize(n - 1);
    sigma.resize(n - 1);

    for (int i = 0; i < n - 1; i++) {
        m_cell[i] = ro0 / 2. * (x[i + 1] - x[i]); // масса ячейки через плотность и длину ячеек
        ro[i] = ro0;
        dsigm[i].resize(3);
        sigma[i].resize(3);
    }


    m_dot.resize(n);

    m_dot[0] = m_cell[0];
    for (int i = 1; i < n - 1; i++) {
        m_dot[i] = m_cell[i - 1] + m_cell[i];
    }
    m_dot[n - 1] = m_cell[n - 2];

    dt = 0;
    t = 0;

    // левый и правый концы закреплены - скорость 0
    is_v0_left = true;
    is_v0_right = true;
    v0_l = 0;
    v0_r = 0;

    // левый и правый концы закреплены - скорость 0
    is_f0_left = false;
    is_f0_right = false;
    f0_l = 0;
    f0_r = 0;

    // установка граничных и нач условий
    set_border_v();
    set_start();

}
void FEM::set_border_v(bool left1, double v0_l1, bool right1, double v0_r1) {
    is_v0_left = left1;
    v0_l = v0_l1;
    is_v0_right = right1;
    v0_r = v0_r1;
}
void FEM::set_border_f(bool left1, double f0_l1, bool right1, double f0_r1) {
    is_f0_left = left1;
    f0_l = f0_l1;
    is_f0_right = right1;
    f0_r = f0_r1;
}
void FEM::set_print_parametres(bool printAnalitcalSolution) {

    printAS = printAnalitcalSolution;

    ffu.open("output/uOut.csv", ios::out);
    ffs.open("output/sOut.csv", ios::out);
    ffu << "Time;\"U(l)\";\"U(l/2)\";";
    if (printAS) {
        ffu << "\"U(l)A\";\"U(l/2)A\"";
    }
    ffu << endl;

    ffs << "Time;\"S(0l)\";\"S(l/2)\";";
    if (printAS) {
        ffs << "\"S(0l)An\";\"S(l/2)An\"";
    }
    ffs << endl;

    cout << "t = ";
}

