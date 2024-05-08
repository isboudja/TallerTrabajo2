//
// Created by isboudja on 08/05/2024.
//


#include <iostream>
#include "Matrix.h"
#include "SAT_Const.h"
#include "Mjday.h"
#include "Position.h"
#include "globals.h"

struct AuxParam {
    double Mjd_UTC;
    int n;
    int m;
    int planets;
    int sun;
    int moon;
};

int main() {

    // Carga de archivos
    // Código para cargar archivos GGM03S.txt, eop19620101.txt, y GEOS3.txt
    AuxParam auxParam;
    int nobs = 46;

    FILE *fid = fopen("../texts/GEOS3.txt","r");
    Matrix obs(nobs,4),Rs(1,3);
    double Y;
    double M;
    double d;
    double H;
    double min;
    double s;
    double az;
    double el;
    double Dist;
    if(fid==nullptr){
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n=0;n<=46;n++){
        fscanf(fid, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",&Y,&M,&d,&H,&min,&s,&az,&el,%Dist);
        obs(n,1) = Mjday(Y,M,d,H,m,s);
        obs(n,2) = Constants::Rad*az;
        obs(n,3) = Constants::Rad*el;
        obs(n,4) = 1e3*Dist;

    }

    fclose(fid);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224 * Constants::Rad; // [rad]
    double sigma_el = 0.0139 * Constants::Rad; // [rad]


    double lat = Constants::Rad * 21.5748;   // [rad]
    double lon = Constants::Rad * (-158.2706); // [rad]
    double alt = 300.20;                // [m]


    Position(lon, lat, alt, Rs);


    double P[6][6] = {{0}};
    double LT_matrix[3][3] = {{0}};



    double Y0_apr[6] = {0};

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = obs(9,1);
    auxParam.sun = 1;
    auxParam.moon = 1;
    auxParam.planets = 1;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.Mjd_UTC = Mjd_UTC;

    for (int i = 0; i < nobs; ++i) {

    }


    double Y_true[6] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};


    double error_pos[3] = {Y0[0] - Y_true[0], Y0[1] - Y_true[1], Y0[2] - Y_true[2]};
    double error_vel[3] = {Y0[3] - Y_true[3], Y0[4] - Y_true[4], Y0[5] - Y_true[5]};


    std::cout << "Error of Position Estimation" << std::endl;
    std::cout << "dX: " << error_pos[0] << " [m]" << std::endl;
    std::cout << "dY: " << error_pos[1] << " [m]" << std::endl;
    std::cout << "dZ: " << error_pos[2] << " [m]" << std::endl;
    std::cout << "Error of Velocity Estimation" << std::endl;
    std::cout << "dVx: " << error_vel[0] << " [m/s]" << std::endl;
    std::cout << "dVy: " << error_vel[1] << " [m/s]" << std::endl;
    std::cout << "dVz: " << error_vel[2] << " [m/s]" << std::endl;

    return 0;
}