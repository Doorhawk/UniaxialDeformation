﻿#include <iostream>
#include <fstream>
#include "FiniteElementMethod.h"


int main()
{
    fstream fin("input/input.txt");

    int n=1;
    double Tmax = 1;
    double l = 0, ro0 = 0, K = 0, G = 0;
    bool leftF = false, rightF = false, leftV = false, rightV = false;
    double f0l = 0, f0r = 0, v0l = 0, v0r = 0;
    bool printAnalit = false;

    string s;
    fin >> s >> n;
    fin >> s >> Tmax;
    fin >> s >> l;
    fin >> s >> ro0;
    fin >> s >> K;
    fin >> s >> G;
    fin >> s >> leftF;
    fin >> s >> rightF;
    fin >> s >> leftV;
    fin >> s >> rightV;
    fin >> s >> f0l;
    fin >> s >> f0r;
    fin >> s >> v0l;
    fin >> s >> v0r;
    fin >> s >> printAnalit;



    if (n <=1 || Tmax <= 0 || l <= 0 || ro0 <= 0 || K <= 0 || G <= 0) {
        cout << "Error - parameter must be > 0 and n > 1";
        return -1;
    }
    if (leftF && leftV) {
        cout << "Error - only one condition is allowed - leftF or leftV";
        return -1;
    }
    if (rightF && rightV) {
        cout << "Error - only one condition is allowed - rightF or rightV";
        return -1;
    }
    if(printAnalit){
        if (!(rightF && leftV && !v0l)) {
            cout << "Error - you are trying printAnalit, it is possible only if rightF = 1 && leftV = 1 and v0l = 0, f0r!=0";
            return -1;
        }
    }

    FEM fem(n,l,ro0,K,G);
    fem.set_border_f(leftF, f0l, rightF, f0r);
    fem.set_border_v(leftV, v0l, rightV, v0r);
    fem.set_print_parametres(printAnalit);

    fem.solve(Tmax);

    fin.close();

    cout <<endl<< "Press Enter"<<endl;
    cin.get();

    return 0;
}


