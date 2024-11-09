#include <iostream>
#include <vector>
#include <cmath>

//Including vortex simulation header
#include <ranges>

#include "PointVortexSim.h"


int main() {
    using namespace std;
    //Intitial conditions
    //list of point vortices Initialization
    vector<PtVort> VortexList;

    //Creating point vortices
    PtVort vort1(-2,0,5);
    VortexList.push_back(vort1);
    PtVort vort2(2,0,-5);
    VortexList.push_back(vort2);

    cout<<"Number of Vortices: "<<VortexList.size()<<endl;

    cout<<" "<< VortexList[1].x[0]<<endl;


    //Vort. Dependent Initial Conditions
    int VortexNum = VortexList.size();  //Number of Vortices
    const double DeltaTime = 0.001;    // Time step size (time resolution)
    const int NumSteps = 5000;      // Number of time steps


    //

    //Runge-Kutta IVP solution

    //Time step loop

    for(int t = 0; t < NumSteps; t++) {

        //loop for each K value
        //K1
        for (int i = 0; i < VortexNum; i++) {
            for (int j = 0; j < VortexNum; j++) {
                if (i != j) {
                    double xi = VortexList[i].x[t];
                    double yi = VortexList[i].y[t];
                    double xj = VortexList[j].x[t];
                    double yj = VortexList[j].y[t];
                    double Sj = VortexList[j].S;

                    VortexList[i].k1[0] = dx(xi, xj, yi, yj, Sj) + VortexList[i].k1[0];
                    VortexList[i].k1[1] = dy(xi, xj, yi, yj, Sj) + VortexList[i].k1[1];
                    cout << VortexList[i].k1[0] << "," << VortexList[i].k1[1]<< endl;

                }
            }
        };

        //K2
        for (int i = 0; i < VortexList.size(); i++) {
            for (int j = 0; j < VortexList.size(); j++) {
                if (i != j) {
                    // Positions for point
                    double xi = VortexList[i].x[t];
                    double yi = VortexList[i].y[t];
                    // position for influencer
                    double xj = VortexList[j].x[t];
                    double yj = VortexList[j].y[t];
                    // Influencer Circulation Strength
                    double Sj = VortexList[j].S;

                    // K1 values for each point
                    double k1_xi = VortexList[i].k1[0];
                    double k1_xj = VortexList[j].k1[0];
                    double k1_yi = VortexList[i].k1[1];
                    double k1_yj = VortexList[j].k1[1];

                    VortexList[i].k2[0] = dx(xi + k1_xi*DeltaTime/2,
                        xj + k1_xj*DeltaTime/2,
                        yi + k1_yi*DeltaTime/2,
                        yj + k1_yj*DeltaTime/2,
                        Sj) + VortexList[i].k2[0];

                    VortexList[i].k2[1] = dy(xi + k1_xi*DeltaTime/2,
                        xj + k1_xj*DeltaTime/2,
                        yi + k1_yi*DeltaTime/2,
                        yj + k1_yj*DeltaTime/2,
                        Sj) + VortexList[i].k2[1];
                }
            }
        };

        //K3
        for (int i = 0; i < VortexList.size(); i++) {
            for (int j = 0; j < VortexList.size(); j++) {
                if (i != j) {
                    double xi = VortexList[i].x[t];
                    double yi = VortexList[i].y[t];
                    double xj = VortexList[j].x[t];
                    double yj = VortexList[j].y[t];
                    double Sj = VortexList[j].S;

                    //k2 Values for each point
                    double k2_xi = VortexList[i].k2[0];
                    double k2_xj = VortexList[j].k2[1];
                    double k2_yi = VortexList[i].k2[0];
                    double k2_yj = VortexList[j].k2[1];

                    VortexList[i].k3[0] = dx(xi + k2_xi*DeltaTime/2,
                        xj + k2_xj*DeltaTime/2,
                        yi + k2_yi*DeltaTime/2,
                        yj+ k2_yj*DeltaTime/2,
                        Sj) + VortexList[i].k3[0];

                    VortexList[i].k3[1] = dy(xi + k2_xi*DeltaTime/2,
                        xj + k2_xj*DeltaTime/2,
                        yi + k2_yi*DeltaTime/2,
                        yj+ k2_yj*DeltaTime/2,
                        Sj) + VortexList[i].k3[1];
                }
            }
        };

        //K4
        for (int i = 0; i < VortexList.size(); i++) {
            for (int j = 0; j < VortexList.size(); j++) {
                if (i != j) {
                    double xi = VortexList[i].x[t];
                    double yi = VortexList[i].y[t];
                    double xj = VortexList[j].x[t];
                    double yj = VortexList[j].y[t];
                    double Sj = VortexList[j].S;

                    double k3_xi = VortexList[i].k3[0];
                    double k3_xj = VortexList[j].k3[0];
                    double k3_yi = VortexList[i].k3[1];
                    double k3_yj = VortexList[j].k3[1];

                    VortexList[i].k4[0] = dx(xi + k3_xi*DeltaTime,
                        xj + k3_xj*DeltaTime,
                        yi + k3_yi*DeltaTime,
                        yj + k3_yj*DeltaTime ,
                        Sj) + VortexList[i].k4[0];
                }
            }
        };

        for (int i = 0; i<VortexNum; i++) {
            double xi = VortexList[i].x[t] +
                DeltaTime/6*(VortexList[i].k1[0] + 2*VortexList[i].k2[0] +
                    2*VortexList[i].k3[0] + VortexList[i].k4[0]);
            double yi = VortexList[i].y[t] +
                DeltaTime/6*(VortexList[i].k1[1] + 2*VortexList[i].k2[1] +
                    2*VortexList[i].k3[1] + VortexList[i].k4[1])                ;
            VortexList[i].PosUpdate(xi, yi, t);
        };
    };



    return 0;
};

