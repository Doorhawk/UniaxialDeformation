
#include "AnalyticalSolution.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <iostream>

using namespace std;

class FEM {

private:

    double K;        //  K      - ������ ��������� ������ 
    double G;        //  G      - ������ ������
    double E;        //  E      - ������ ����
    double Nu;       //  Nu     - ��� ��������
    double c;        //  c      - �������� ����� � �����
    int n;           //  n      - ���������� �����, (n-1) - ���������� �����
    double l;        //  l      - ������ ������
    double ro0;      //  ro0    - ��������� ���������
    double t;        //  t      - �����
    double dt;       //  dt     - ��� �� ������� ���������� ����� ti+1 - ti
    double dt2;      //  dt2    - ������������ ��� �� ������� ���������� ����� ti+1/2 - ti-1/2
    vector<double> m_cell;  //  m_cell - ����� ������; 
    vector<double> m_dot;   //  m_dot  - ������� �����
    vector<double> ro;      //  ro     - ��������� ������ ����������
    vector<double> x0;      //  x0     - ��������� ���������� �����
    vector<double> x;       //  x      - ���������� ����� � ����� ����� �� �������
    vector<double> x2;      //  x2     - ��������� ����� � ����������(�������������) ����� �� �������
    vector<double> v;       //  v      - �������� ����
    vector<double> u;       //  u      - ����������� ����
    vector<double> f;       //  f      - ���� � �����
    vector<double> psigm;   //  psigm  - ������� ������ ����������
    vector<vector<double>> dsigm;  //  dsigm  - �������� ������� ����������
    vector<vector<double>> sigma;  //  sigma  - ������ ������ ����������

    bool is_v0_left, is_v0_right;  //  ���� �� ������� ������� �� �����/������ ������ �� ��������
    double v0_l, v0_r;             //  ���� ���� �� ����� �� ��� ������ �������� ��������

    bool is_f0_left, is_f0_right;  //  ���� �� ������� ������� �� �����/������ ������ �� ����
    double f0_l, f0_r;             //  ���� ���� �� ����� �� ��� ������ ����

    bool printAS; // �������� �� ������������� �������

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