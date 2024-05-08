//
// Created by isboudja on 08/05/2024.
//

#include "Accel.h"


/*Matrix Accel(double x,double Y){

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;


    // Cálculo de los parámetros de tiempo
    IERS::GetParams(eopdata, Mjd_UTC + x / 86400, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    double UT1_TAI = UT1_UTC - TAI_UTC;
    double UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = Mjd_UTC + x / 86400 + UT1_UTC / 86400;
    double Mjd_TT = Mjd_UTC + x / 86400 + TT_UTC / 86400;

    // Matrices de transformación
    Matrix3d P = PrecMatrix(const.MJD_J2000, Mjd_TT);
    Matrix3d N = NutMatrix(Mjd_TT);
    Matrix3d T = N * P;
    Matrix3d E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    // Tiempo barycentrico dinamico modificado
    double MJD_TDB = Mjday_TDB(Mjd_TT);

    // Posiciones de los cuerpos celestes
    Vector3d r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun;
    JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

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





}*/