//
// Created by whatl on 11/8/2024.
//
#ifndef POINTVORTEXSIM_H
#define POINTVORTEXSIM_H
#include <iostream>
#include <vector>
#include <cmath>

//Namespace declaration
using namespace std;

//counting number of objects for "PtVort" Class
static int VortexNumb = 0;

//"PtVort" Class, for tracking Vortices
class PtVort {
    public:
    //position and time vectors, for all previous positions
    vector<double> x {};    //X position vector
    vector<double> y {};    //Y position vector
    vector<int> t_ind {};   //time index vector, integer qty

    //Scalar Quantities
    double S;
    int VortexIndex;

    //Runge-Kutta K values
    vector<double> k1 {0,0};
    vector<double> k2 {0,0};
    vector<double> k3 {0,0};
    vector<double> k4 {0,0};

    PtVort(double Initial_X, double Initial_Y, double Circulation_Strength) {

        x.push_back( Initial_X );
        y.push_back( Initial_Y );
        t_ind.push_back(0);
        S = Circulation_Strength;
        VortexIndex = VortexNumb;

        cout << "Initial X: " << x[0] << " | Initial Y: " << y[0] << " | Circulation Strength: " << S << " | Vortex Index: "<< VortexIndex <<endl;

        ++VortexNumb;
    };

    void PosUpdate(double New_X, double New_Y, int Time_Step) {
        x.push_back(New_X);   //Updating X position
        y.push_back(New_Y);   //Updating Y position
        t_ind.push_back(Time_Step);   //Updating time step

        cout<< "New X:" << New_X << " | New Y:" << New_Y << " | time step:" << Time_Step <<endl;
    };
};

inline double dx(double xi, double xj, double yi, double yj, double Sj) {
    double del_x = (xi - xj);
    double del_y = (yi - yj);
    return -1*Sj*(del_y)/(2*numbers::pi*(pow(del_x,2) + pow(del_y,2)));
};

inline double dy(double xi, double xj, double yi, double yj, double Sj) {
    double del_x = (xi - xj);
    double del_y = (yi - yj);
    return Sj*(del_x)/(2*numbers::pi*(pow(del_x,2) + pow(del_y,2)));
};




#endif //POINTVORTEXSIM_H

