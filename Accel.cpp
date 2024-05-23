//
// Created by isboudja on 08/05/2024.
//

#include <cstdio>
#include "Accel.h"
#include "Matrix.h"
#include "SAT_Const.h"
#include "timediff.h"
#include "Mjday_TDB.h"
#include "AuxParam.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "IERS.h"
#include "globals.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"

Matrix Accel(double x,Matrix &Y){

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    AuxParam AuxParam;
    Constants c;
    AuxParam.Mjd_UTC = 4.974611635416653e+04;
    Matrix eopdata(13, 21413);
    FILE *fid = fopen("../texts/eop19620101.txt", "r");

    if (fid == nullptr) {
        printf("error globals");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= 21413; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &eopdata(1, i),
               &eopdata(2, i), &eopdata(3, i),
               &eopdata(4, i), &eopdata(5, i), &eopdata(6, i),
               &eopdata(7, i), &eopdata(8, i), &eopdata(9, i),
               &eopdata(10, i), &eopdata(11, i), &eopdata(12, i),
               &eopdata(13, i));
    }


    fclose(fid);
    // Cálculo de los parámetros de tiempo
    Matrix result(1,9);
    result = IERS(eopdata, AuxParam.Mjd_UTC + x/86400,'n');
    x_pole = result(1,1);
    y_pole = result(1,2);
    UT1_UTC = result(1,3);
    LOD = result(1,4);
    dpsi = result(1,5);
    deps=  result(1,6);
    dx_pole=result(1,7) ;
    dy_pole=result(1,8) ;
    TAI_UTC = result(1,9);
    double UT1_TAI = UT1_UTC - TAI_UTC;
    double UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400 + UT1_UTC / 86400;
    double Mjd_TT = AuxParam.Mjd_UTC + x / 86400 + TT_UTC / 86400;

    // Matrices de transformación
    Matrix P = PrecMatrix(c.MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    // Tiempo barycentrico dinamico modificado
    double MJD_TDB = Mjday_TDB(Mjd_TT);

    Matrix r_Mercury(1,3), r_Venus(1,3), r_Earth(1,3), r_Mars(1,3), r_Jupiter(1,3), r_Saturn(1,3), r_Uranus(1,3), r_Neptune(1,3), r_Pluto(1,3), r_Moon(1,3), r_Sun(1,3);
    Matrix result2(1,33);
    result2 = JPL_Eph_DE430(MJD_TDB);
    for(int j=1;j<=3;j++){
        r_Mercury(1,j) = result2(1,j*1);
        r_Venus(1,j) = result2(1,j+3);
        r_Earth(1,j) = result2(1,j+6);
        r_Mars(1,j) = result2(1,j+9);
        r_Jupiter(1,j) = result2(1,j+12);
        r_Saturn(1,j) = result2(1,j+15);
        r_Uranus(1,j) = result2(1,j+18);
        r_Neptune(1,j) = result2(1,j+21);
        r_Pluto(1,j) = result2(1,j+24);
        r_Moon(1,j) = result2(1,j+27);
        r_Sun(1,j) = result2(1,j+30);

    }
    // Aceleración debida al campo gravitacional armónico
    Matrix a(1,1);
    Matrix r(1,3);
    Matrix a2(1,3);
    for(int i=1;i<=3;i++){
        r(1,i) = Y(1,i);
        a2(1,i) = 1;
    }
    a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    a2*a;
    a2.print();
    // Perturbaciones luni-solares
    if (AuxParam.sun) {
        a = a + AccelPointMass(r, r_Sun, c.GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r, r_Moon, c.GM_Moon);
    }

    // Perturbaciones planetarias
    if (AuxParam.planets) {
        a = a + AccelPointMass(r, r_Mercury, c.GM_Mercury);
        a = a + AccelPointMass(r, r_Venus, c.GM_Venus);
        a = a + AccelPointMass(r, r_Mars, c.GM_Mars);
        a = a + AccelPointMass(r, r_Jupiter, c.GM_Jupiter);
        a = a + AccelPointMass(r, r_Saturn, c.GM_Saturn);
        a = a + AccelPointMass(r, r_Uranus, c.GM_Uranus);
        a = a + AccelPointMass(r, r_Neptune, c.GM_Neptune);
        a = a + AccelPointMass(r, r_Pluto, c.GM_Pluto);
    }

    Matrix dY(1,4);
    for(int i=1;i<=3;i++){
    dY(1,i) = Y(1,i+3);
    }
    dY(1,4) = a(1,1);
    dY.print();
    return dY;

}





