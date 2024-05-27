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
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "JPL_Eph_DE430.h"
#include "globals.h"

Matrix Accel(double x,Matrix &Y){

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    AuxParam AuxParam;
    Constants c;
    AuxParam.Mjd_UTC = 4.974611635416653e+04;


    // Cálculo de los parámetros de tiempo

    Matrix result(1,9);
    result = IERS(*globals::eopdata, AuxParam.Mjd_UTC + x/86400,'n');
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
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;// Tiempo barycentrico dinamico modificado
    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury2(1,3);
    Matrix r_Venus2(1,3);
    Matrix r_Earth2(1,3) ;
    Matrix r_Mars2(1,3);
    Matrix r_Jupiter2(1,3);
    Matrix r_Saturn2(1,3);
    Matrix r_Uranus2(1,3);
    Matrix r_Neptune2(1,3);
    Matrix r_Pluto2(1,3);
    Matrix r_Moon2(1,3);
    Matrix r_Sun2(1,3);
    Matrix  res = JPL_Eph_DE430(MJD_TDB);
    for(int j=1;j<=3;j++){
        r_Mercury2(1,j) = res(1,j*1);
        r_Venus2(1,j) = res(1,j+3);
        r_Earth2(1,j) = res(1,j+6);
        r_Mars2(1,j) = res(1,j+9);
        r_Jupiter2(1,j) = res(1,j+12);
        r_Saturn2(1,j) = res(1,j+15);
        r_Uranus2(1,j) = res(1,j+18);
        r_Neptune2(1,j) = res(1,j+21);
        r_Pluto2(1,j) = res(1,j+24);
        r_Moon2(1,j) = res(1,j+27);
        r_Sun2(1,j) = res(1,j+30);

    }
    // Aceleración debida al campo gravitacional armónico

    Matrix r(3,1);
    Matrix a2(1,3);
    for(int i=1;i<=3;i++){
        r(i,1) = Y(i,1);
        a2(1,i) = 1;
    }
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;
    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    Matrix r2 = r.transpose();


    // Perturbaciones luni-solares
    if (AuxParam.sun) {
        a = a + AccelPointMass(r2, r_Sun2, c.GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r2, r_Moon2, c.GM_Moon);
    }

    // Perturbaciones planetarias
    if (AuxParam.planets) {
        a = a + AccelPointMass(r2, r_Mercury2, c.GM_Mercury);
        a = a + AccelPointMass(r2, r_Venus2, c.GM_Venus);
        a = a + AccelPointMass(r2, r_Mars2, c.GM_Mars);
        a = a + AccelPointMass(r2, r_Jupiter2, c.GM_Jupiter);
        a = a + AccelPointMass(r2, r_Saturn2, c.GM_Saturn);
        a = a + AccelPointMass(r2, r_Uranus2, c.GM_Uranus);
        a = a + AccelPointMass(r2, r_Neptune2, c.GM_Neptune);
        a = a + AccelPointMass(r2, r_Pluto2, c.GM_Pluto);
    }

    Matrix dY(1,6);
    for(int i=1;i<=3;i++){
    dY(1,i) = Y(i+3,1);
    }
    dY(1,4) = a(1,1);
    dY(1,5) = a(1,2);
    dY(1,6) = a(1,3);
    return dY;

}





