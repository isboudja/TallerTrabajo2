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

Matrix Accel(double x,double Y){

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    AuxParam AuxParam;
    Constants c;

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
    // Cálculo de los parámetros de tiempo
    Matrix result(1,9);
    result = IERS(eopdata, AuxParam.Mjd_UTC + x / 86400,'n');
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

    // Posiciones de los cuerpos celestes
    double r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun;
    JPL_Eph_DE430(MJD_TDB);

    // Aceleración debida al campo gravitacional armónico
    Vector3d a = AccelHarmonic(Y.head(3), E, AuxParam.n, AuxParam.m);

    // Perturbaciones luni-solares
    if (AuxParam.sun) {
        a += AccelPointMass(Y.head(3), r_Sun, const.GM_Sun);
    }
    if (AuxParam.moon) {
        a += AccelPointMass(Y.head(3), r_Moon, const.GM_Moon);
    }

    // Perturbaciones planetarias
    if (AuxParam.planets) {
        a += AccelPointMass(Y.head(3), r_Mercury, const.GM_Mercury);
        a += AccelPointMass(Y.head(3), r_Venus, const.GM_Venus);
        a += AccelPointMass(Y.head(3), r_Mars, const.GM_Mars);
        a += AccelPointMass(Y.head(3), r_Jupiter, const.GM_Jupiter);
        a += AccelPointMass(Y.head(3), r_Saturn, const.GM_Saturn);
        a += AccelPointMass(Y.head(3), r_Uranus, const.GM_Uranus);
        a += AccelPointMass(Y.head(3), r_Neptune, const.GM_Neptune);
        a += AccelPointMass(Y.head(3), r_Pluto, const.GM_Pluto);
    }

    Vector6d dY;
    dY << Y.tail(3), a;
    return dY;
}





}