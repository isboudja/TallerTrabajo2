//
// Created by isboudja on 08/05/2024.
//


#include <iostream>
#include "Matrix.h"
#include "SAT_Const.h"
#include "Mjday.h"
#include "Position.h"
#include "globals.h"
#include "IERS.h"
#include "timediff.h"
#include "gmst.h"
#include "R_z.h"
#include "AzElPa.h"
#include "AuxParam.h"


int main() {

    // Carga de archivos
    // Código para cargar archivos GGM03S.txt, eop19620101.txt, y GEOS3.txt

    Matrix eopdata(13, 21413);
    FILE *fid = fopen("../texts/eop19620101.txt", "r");

    if (fid == nullptr) {
        printf("error globals");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= 21413; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &(*globals::matrix)(1, i),
               &(*globals::matrix)(2, i), &(*globals::matrix)(3, i),
               &(*globals::matrix)(4, i), &(*globals::matrix)(5, i), &(*globals::matrix)(6, i),
               &(*globals::matrix)(7, i), &(*globals::matrix)(8, i), &(*globals::matrix)(9, i),
               &(*globals::matrix)(10, i), &(*globals::matrix)(11, i), &(*globals::matrix)(12, i),
               &(*globals::matrix)(13, i));
    }


    fclose(fid);

    AuxParam auxParam;
    int nobs = 46;

    FILE *fid2 = fopen("../texts/GEOS3.txt", "r");
    Matrix obs(nobs, 4), Rs(1, 3);
    double y;
    double M;
    double d;
    double H;
    double min;
    double s;
    double az;
    double el;
    double Dist;
    if (fid == nullptr) {
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n = 0; n <= 46; n++) {
        fscanf(fid2, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &y, &M, &d, &H, &min, &s, &az, &el, &Dist);
        obs(n, 1) = Mjday(y, M, d, H, m, s);
        obs(n, 2) = Constants::Rad * az;
        obs(n, 3) = Constants::Rad * el;
        obs(n, 4) = 1e3 * Dist;

    }

    fclose(fid2);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224 * Constants::Rad; // [rad]
    double sigma_el = 0.0139 * Constants::Rad; // [rad]


    double lat = Constants::Rad * 21.5748;   // [rad]
    double lon = Constants::Rad * (-158.2706); // [rad]
    double alt = 300.20;                // [m]


    Position(lon, lat, alt, Rs);
    int Mjd1 = obs(1, 1);
    int Mjd2 = obs(9, 1);
    int Mjd3 = obs(18, 1);

    double LT_matrix[3][3] = {{0}};


    Matrix Y0_apr(6, 6);

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = obs(9, 1);
    auxParam.sun = 1;
    auxParam.moon = 1;
    auxParam.planets = 1;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.Mjd_UTC = Mjd_UTC;


    double n_eqn = 6;
    Matrix Y(6,1);
    //Y = DEInteg(@Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    Matrix P(6, 6);

    for (int i = 1; i < 4; i++) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i < 7; i++) {
        P(i, i) = 1e3;
    }


    //LT = LTC(lon,lat);
    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    int t = 0;

    for (int i = 1; i <= nobs; i++) {
        int t_old = t;
        Matrix Y_old(6, 1) = Y;

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;
        Matrix Resulteop(1, 9);
        Resulteop = IERS(eopdata, Mjd_UTC, 'l');
        double x_pole = Resulteop(1, 1);
        double y_pole = Resulteop(1, 2);
        double UT1_UTC = Resulteop(1, 3);
        double LOD = Resulteop(1, 4);
        double dpsi = Resulteop(1, 5);
        double deps = Resulteop(1, 6);
        double dx_pole = Resulteop(1, 7);
        double dy_pole = Resulteop(1, 8);
        double TAI_UTC = Resulteop(1, 9);

        double UT1_TAI;
        double UTC_GPS;
        double UT1_GPS;
        double TT_UTC;
        double GPS_UTC;
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        auxParam.Mjd_UTC = Mjd_UTC;
        auxParam.Mjd_TT = Mjd_TT;


    for (int ii = 1; ii <= 6; ii++) {
        yPhi(ii, 1) = Y_old(ii, 1);
        for (int j = 1; j <= 6; j++){
            if (ii == j) {
                yPhi(6 * j + ii,1) = 1;
            }else{
                yPhi(6 * j + ii,1) = 0;
    }
}
}
    // yPhi = DEInteg (@VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

    for (int j = 1; j <= 6; j++) {
        for(int nnn=1; nnn <Phi.fil;nnn++) {
            Phi(nnn,j) = yPhi(6 * j + nnn,1);
        }
    }

    //Y = DEInteg (@Accel,0,t-t_old,1e-13,1e-6,6,Y_old);


        double theta = gmst(Mjd_UT1);
        Matrix U(3,3);
        U = R_z(theta);
        Matrix r(1,3);
        r(1,1) = Y(1,1);
        r(1,2) = Y(1,2);
        r(1,3) = Y(1,3);
        Matrix LT(3,3);
        Matrix s(1,3);
        s = LT*(U*r-Rs);


        //P = TimeUpdate(P, Phi);

        Matrix AZEL(1,8);
        AZEL = AzElPa(s);
        dAdY = [dAds*LT*U,zeros(1,3)];


        [K, Y, P] = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );


        r = Y(1:3);
        s = LT*(U*r-Rs);
        [Azim, Elev, dAds, dEds] = AzElPa(s);
        dEdY = [dEds*LT*U,zeros(1,3)];


        [K, Y, P] = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );


        r = Y(1:3);
        s = LT*(U*r-Rs);
        Dist = norm(s); dDds = (s/Dist)';
        dDdY = [dDds*LT*U,zeros(1,3)];

        % Measurement update
        [K, Y, P] = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
        end

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