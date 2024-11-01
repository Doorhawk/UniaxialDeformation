
#include "AnalyticalSolution.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <iostream>

using namespace std;

class FEM {

private:

    double K;        //  K      - модуль обьемного сжати€ 
    double G;        //  G      - модуль сдвига
    double E;        //  E      - модуль юнга
    double Nu;       //  Nu     - кэф пуассона
    double c;        //  c      - скорость звука в среде
    int n;           //  n      - количество узлов, (n-1) - количество €чеек
    double l;        //  l      - длинна пр€мой
    double ro0;      //  ro0    - начальна€ плотность
    double t;        //  t      - врем€
    double dt;       //  dt     - шаг по времени рассто€ние между ti+1 - ti
    double dt2;      //  dt2    - половинчатый шаг по времени рассто€ние между ti+1/2 - ti-1/2
    vector<double> m_cell;  //  m_cell - масса €чейки; 
    vector<double> m_dot;   //  m_dot  - узлова€ амсса
    vector<double> ro;      //  ro     - плотность €чейки актуальна€
    vector<double> x0;      //  x0     - начальные координаты узлов
    vector<double> x;       //  x      - координаты узлов с целым шагом по времени
    vector<double> x2;      //  x2     - кординаты узлов с половинным(промежуточным) шагом по времени
    vector<double> v;       //  v      - скорость узла
    vector<double> u;       //  u      - перемешешие узла
    vector<double> f;       //  f      - силы в узлах
    vector<double> psigm;   //  psigm  - шаровой тензор напр€жений
    vector<vector<double>> dsigm;  //  dsigm  - девиатор тензора напр€жений
    vector<vector<double>> sigma;  //  sigma  - полный тензор напр€жений

    bool is_v0_left, is_v0_right;  //  есть ли краевые услови€ на левом/правом конеце по — ќ–ќ—“»
    double v0_l, v0_r;             //  если есть то кака€ на нем задана скорость скорость

    bool is_f0_left, is_f0_right;  //  есть ли краевые услови€ на левом/правом конеце по —»Ћ≈
    double f0_l, f0_r;             //  если есть то кака€ на нем задана сила

    bool printAS; // печатать ли аналитическое решение

    AnalitcalSolution ana;
    fstream ffu;
    fstream ffs;

    void set_border_v();
    void set_border_f();
    void set_start();

    void calc_dti_ro();
    void calc_x_u();
    void calc_eps_sigm();
    void calc_forse();
    void calc_v();
    void correct_v();

    void printSelectedInformation();
    void step();
public:
    FEM(int _n, double _l, double _ro0, double _K, double _G);

    void set_border_v(bool left1, double v0_l1, bool right1, double v0_r1);
    void set_border_f(bool left1, double f0_l1, bool right1, double f0_r1);
    void set_print_parametres(bool printAnalitcalSolution);

    void solve(double Tmax);

};