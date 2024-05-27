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
#include "Accel.h"
#include "DEInteg.h"
#include "VarEqn.h"
#include "LTC.h"
#include "TimeUpdate.h"
#include "MeasUpdate.h"
#include "globals.h"
/*
int main() {




    globals::eop1962();
    globals::DE430();
    globals::GGM();

    AuxParam auxParam;
    int nobs = 46;
    globals::GEOS3(nobs);
    double y;
    double M;
    double d;
    double H;
    double min;
    double s;
    double az;
    double el;
    double Dist;
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    int n;
    int m;
    int ii;
    int i;
    int j;
    int l;

    Matrix Rs(nobs, 4);


    double sigma_range = 92.5;
    double sigma_az = 0.0224 * Constants::Rad;
    double sigma_el = 0.0139 * Constants::Rad;


    double lat = Constants::Rad * 21.5748;
    double lon = Constants::Rad * (-158.2706);
    double alt = 300.20;


    Position(lon, lat, alt, Rs);
    int Mjd1 = (*globals::obs)(1, 1);
    int Mjd2 = (*globals::obs)(9, 1);
    int Mjd3 = (*globals::obs)(18, 1);


    Matrix LT_matrix(3,3);

    Matrix Y0_apr(6, 1);
    Y0_apr(1,1) =  6221397.62857869;
    Y0_apr(2,1) = 2867713.77965738 ;
    Y0_apr(3,1) = 3006155.98509949;
    Y0_apr(4,1) = 4645.04725161807;
    Y0_apr(5,1) = -2752.21591588205;
    Y0_apr(6,1) =  -7507.99940987033;


    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = (*globals::obs)(9, 1);
    auxParam.sun = 1;
    auxParam.moon = 1;
    auxParam.planets = 1;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.Mjd_UTC = Mjd_UTC;


    double n_eqn = 6;
    Matrix Y(6,1);
    Y0_apr.print();

    DEInteg(Accel,0,-((*globals::obs)(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
    Y = Y0_apr;
    Matrix P(6, 6);

    for (i = 1; i < 4; i++) {
        P(i, i) = 1e8;
    }
    for (i = 4; i < 7; i++) {
        P(i, i) = 1e3;
    }


    Matrix LT = LTC(lon,lat);
    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    int t = 0;
    int t_old;
    Matrix Y_old(6,1);
    Matrix Resulteop(1, 9);
    Matrix AZEL(1,8);
    Matrix Azim(1,2);
    Matrix Elev(1,2);
    Matrix dAds(1,2);
    Matrix dEds(1,2);
    Matrix dAdY2(1,3);
    Matrix cer(1,3);
    Matrix dAdY(1,6);
    Matrix dis(1,1);
    Matrix K(6,1);
    Matrix dEdY2(1,3);
    Matrix dEdY(1,6);
    Matrix U(3,3);
    Matrix r(1,3);
    Matrix s1(1,3);
    double Mjd_TT;
    double Mjd_UT1;
    double theta;
    Matrix dDds(3,1);
    Matrix DD(1,3);
    Matrix dDdY2(1,3);
    Matrix dDdY(1,6);
    for (i = 1; i <= nobs; i++) {
        t_old = t;
        Y_old = Y;

        Mjd_UTC = (*globals::obs)(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        Resulteop = IERS(*globals::eopdata, Mjd_UTC, 'l');
        x_pole = Resulteop(1, 1);
        y_pole = Resulteop(1, 2);
         UT1_UTC = Resulteop(1, 3);
        LOD = Resulteop(1, 4);
         dpsi = Resulteop(1, 5);
         deps = Resulteop(1, 6);
        dx_pole = Resulteop(1, 7);
         dy_pole = Resulteop(1, 8);
         TAI_UTC = Resulteop(1, 9);


        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

         Mjd_TT = Mjd_UTC + TT_UTC / 86400;

         Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        auxParam.Mjd_UTC = Mjd_UTC;
        auxParam.Mjd_TT = Mjd_TT;


    for ( ii = 1; ii <= 6; ii++) {
        yPhi(ii, 1) = Y_old(ii, 1);
        for (int j = 1; j <= 6; j++){
            if (ii == j) {
                yPhi(6 * j + ii,1) = 1;
            }else{
                yPhi(6 * j + ii,1) = 0;
    }
}
}
        DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

    for (j = 1; j <= 6; j++) {
        for( n=1; n <Phi.fil;n++) {
            Phi(n,j) = yPhi(6 * j + n,1);
        }
    }

       DEInteg(Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
        Y = Y_old;

        theta = gmst(Mjd_UT1);

        U = R_z(theta);

        r(1,1) = Y(1,1);
        r(1,2) = Y(1,2);
        r(1,3) = Y(1,3);

        s1 = LT*(U*r-Rs);


        P = TimeUpdate(P, Phi,0.0);


        AZEL = AzElPa(s1);
        for(i=1;i<=2;i++){
            Azim(1,i) = AZEL(1,i);
            Elev(1,i) = AZEL(1,i+2);
            dAds(1,i) = AZEL(1,i+4);
            dEds(1,i) = AZEL(1,i+6);
        }


         dAdY2 = dAds*LT*U;

        dAdY = Matrix::concat(dAdY2,cer);

        K = MeasUpdate( Y, (*globals::obs)(i,2), Azim, sigma_az, dAdY, P, 6);


        r(1,1) = Y(1,1);
        r(1,2) = Y(1,2);
        r(1,3) = Y(1,3);
        s1 = LT*(U*r-Rs);
        AZEL = AzElPa(s1);
        for(i=1;i<=2;i++){
            Azim(1,i) = AZEL(1,i);
            Elev(1,i) = AZEL(1,i+2);
            dAds(1,i) = AZEL(1,i+4);
            dEds(1,i) = AZEL(1,i+6);
        }

         dEdY2 = dEds*LT*U;
         dEdY = Matrix::concat(dEdY2,cer);

        K = MeasUpdate ( Y, (*globals::obs)(i,3), Elev, sigma_el, dEdY, P, 6 );
        r(1,1) = Y(1,1);
        r(1,2) = Y(1,2);
        r(1,3) = Y(1,3);
        s1 = LT*(U*r-Rs);
        Dist = s1.norm();

        dis(1,1) = Dist;
         dDds = (s1/Dist);
         DD = dDds.transpose();
         dDdY2 = DD*LT*U;
         dDdY = Matrix::concat(dDdY2,cer);
        K= MeasUpdate( Y, (*globals::obs)(i,4), dis, sigma_range, dDdY, P, 6 );
    }

    Resulteop = IERS(*globals::eopdata, (*globals::obs)(46,1), 'l');
     x_pole = Resulteop(1, 1);
    y_pole = Resulteop(1, 2);
     UT1_UTC = Resulteop(1, 3);
     LOD = Resulteop(1, 4);
     dpsi = Resulteop(1, 5);
     deps = Resulteop(1, 6);
     dx_pole = Resulteop(1, 7);
     dy_pole = Resulteop(1, 8);
     TAI_UTC = Resulteop(1, 9);

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

     Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.Mjd_TT = Mjd_TT;

    DEInteg(Accel,0,-((*globals::obs)(46,1)-(*globals::obs)(1,1))*86400.0,1e-13,1e-6,6,Y);
    Matrix Y0 = Y;
    double Y_t[6] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true(6,1,Y_t,6);

    double error_pos[3] = {Y0(1,1) - Y_true(1,1), Y0(2,1) - Y_true(2,1), Y0(3,1) - Y_true(3,1)};
    double error_vel[3] = {Y0(4,1) - Y_true(4,1), Y0(5,1) - Y_true(4,1), Y0(6,1) - Y_true(6,1)};


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
*/