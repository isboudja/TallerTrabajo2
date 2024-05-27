//
// Created by itsfr on 25/05/2024.
//

#include <cstdio>
#include "VarEqn.h"
#include "AuxParam.h"
#include "SAT_Const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "globals.h"

/**
 * @param x El tiempo de entrada en segundos.
 * @param yPhi La matriz de estado y matriz de covarianza asociada en el tiempo de entrada.
 * @return La derivada de la matriz de estado con respecto al tiempo.
 */

Matrix VarEqn(double x,Matrix &yPhi){

AuxParam AuxParam;
Constants c;
AuxParam.Mjd_UTC = 4.974611635416653e+04;

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;
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
    Matrix r(3,1);
    Matrix v(3,1);
    Matrix yPhip(42,1);
    Matrix Phip(6,6);
    Matrix dfdy(6,6);
    int j;
    int i;
    yPhi.print();
    for(j=1;j<=3;j++){
        r(j,1) = yPhi(j,1);
        v(j,1) = yPhi(j+3,1);
    }
    Matrix Phi(6,6);

    for (j=1;j<=6;j++){
        Phi(1,j) = yPhi(6*j+1,1);
        Phi(2,j) = yPhi(6*j+2,1);
        Phi(3,j) = yPhi(6*j+3,1);
        Phi(4,j) = yPhi(6*j+4,1);
        Phi(5,j) = yPhi(6*j+5,1);
        Phi(6,j) = yPhi(6*j+6,1);
    }

    Matrix a = AccelHarmonic( r, E, AuxParam.n, AuxParam.m);
    Matrix G = G_AccelHarmonic( r, E, AuxParam.n, AuxParam.m);


    for (i=1;i<=3;i++){
    for (j=1;j<=3;j++){
    dfdy(i,j) = 0.0;
    dfdy(i+3,j) = G(i,j);
    if ( i==j ){
        dfdy(i,j+3) = 1;
    }else{
        dfdy(i,j+3) = 0;
        }
        dfdy(i+3,j+3) = 0.0;
    }
    }

    Phip = dfdy*Phi;

    for (i=1;i<=3;i++){
        yPhip(i,1)   = v(i,1);
        yPhip(i+3,1) = a(1,i);
    }

    for (i=1;i<=6;i++){
        for (j=1;j<=6;j++){
            yPhip(6*j+i,1) = Phip(i,j);
        }
    }

    return yPhip;

}