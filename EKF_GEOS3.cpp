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


int main() {

    // Carga de archivos
    // Código para cargar archivos GGM03S.txt, eop19620101.txt, y GEOS3.txt



    AuxParam auxParam;
    int nobs = 46;




    FILE *fid3 = fopen("../texts/GGM03S.txt","r");
    Matrix Cnm(181,181);
    Matrix Snm(181,181);
    Matrix temp(1,6);
    double c;
    if(fid3==nullptr){
        printf("error");
        exit(EXIT_FAILURE);
    }

    for (n=0;n<=180;n++){
        for (m=0;m<=n;m++) {
            fscanf(fid3, "%d%d%lf%lf%lf%lf",&c,&c,&temp(1,3),&temp(1,4),&temp(1,5),&temp(1,6));
            Cnm(n+1,m+1) = temp(1,3);
            Snm(n+1,m+1) = temp(1,4);
        }
    }

    fclose(fid3);

    FILE *fid4 = fopen("../texts/DE430Coeff.txt","r");
    Matrix PC(2285,1020);

    if(fid4==nullptr){
        printf("error2");
        exit(EXIT_FAILURE);
    }

    for (n=1;n<=2285;n++){
        for (m=1;m<=1020;m++) {
            fscanf(fid4, "%lf,",&PC(n,m));
        }
    }


    fclose(fid4);

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


    Matrix LT_matrix(3,3);

    Matrix Y0_apr(6, 1);
    Y0_apr(1,1) =  6221397.62857869;
    Y0_apr(2,1) = 2867713.77965738 ;
    Y0_apr(3,1) = 3006155.98509949;
    Y0_apr(4,1) = 4645.04725161807;
    Y0_apr(5,1) = -2752.21591588205;
    Y0_apr(6,1) =  -7507.99940987033;


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
    DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr,eopdata,Snm,Cnm,PC);
    Y = Y0_apr;
    double gd = -(obs(9,1)-Mjd0)*86400.0;
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

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        Resulteop = IERS(eopdata, Mjd_UTC, 'l');
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
        DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi,eopdata,Snm,Cnm,PC);

    for (j = 1; j <= 6; j++) {
        for( n=1; n <Phi.fil;n++) {
            Phi(n,j) = yPhi(6 * j + n,1);
        }
    }

       DEInteg(Accel,0,t-t_old,1e-13,1e-6,6,Y_old,eopdata,Snm,Cnm,PC);
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

        K = MeasUpdate( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6);


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

        K = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
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
        K= MeasUpdate( Y, obs(i,4), dis, sigma_range, dDdY, P, 6 );
    }

    Resulteop = IERS(eopdata, obs(46,1), 'l');
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

    DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y,eopdata,Snm,Cnm,PC);
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