#include <iostream>
#include <cmath>
#include "vector.h"
#include "seebatt.h"
#include "R_x.h"
#include "globals.h"
#include "Matrix.h"
#include "R_y.h"
#include "R_z.h"
#include "AccelPointMass.h"
#include "AzElPa.h"
#include "Cheb3D.h"
#include "SAT_Const.h"
#include "Frac.h"
#include "Geodetic.h"
#include "Legendre.h"
#include "EccAnom.h"
#include "MeanObliquity.h"
#include "Mjday.h"
#include "Mjday_TDB.h"
#include "NutAngles.h"
#include "Position.h"
#include "sign_.h"
#include "timediff.h"
#include "unit.h"
#include "IERS.h"
#include "globals.h"

#define TOL_ 10e-14

using namespace std;
int tests_run = 0;
const double Constants::pi2 = 2 * M_PI;
const double Constants::Rad = M_PI / 180;
const double Constants::Deg = 180 / M_PI;
const double Constants::Arcs = 3600 * 180 / M_PI;
const double Constants::MJD_J2000 = 51544.5;
const double Constants::T_B1950 = -0.500002108;
const double Constants::c_light = 299792458.000000000;
const double Constants::AU = 149597870700.000000;
const double Constants::f_Earth = 1 / 298.257223563;
const double Constants::R_Sun = 696000e3;
const double Constants::R_Moon = 1738e3;
const double Constants::omega_Earth = 15.04106717866910 / 3600 * Constants::Rad;
const double Constants::GM_Earth = 3.986e14;
const double Constants::GM_Sun = 132712440041.939400e9;
const double Constants::GM_Moon = Constants::GM_Earth / 81.30056907419062;
const double Constants::GM_Mercury = 22031.780000e9;
const double Constants::GM_Venus = 324858.592000e9;
const double Constants::GM_Mars = 42828.375214e9;
const double Constants::GM_Jupiter = 126712764.800000e9;
const double Constants::GM_Saturn = 37940585.200000e9;
const double Constants::GM_Uranus = 5794548.600000e9;
const double Constants::GM_Neptune = 6836527.100580e9;
const double Constants::GM_Pluto = 977.0000000000009e9;
#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int EqMatrix(Matrix &m,Matrix &n,double tol)
{
    int mf = m.fils();
    int nf = n.fils();
    int mc = m.cols();
    int nc = n.cols();


    if(mf != nf || nc != mc){
        printf("false %lf %lf %lf %lf\n",mf, nf, nc, mc);
        return false;
    }else {
        for (int i = 1; i <= mf; i++) {
            for (int j = 1; j <= nc; j++)
                if (fabs(m(i, j) - n(i, j)) > tol) {
                    return false;
                }

        }
    }

    return true;
}

int RX()
{
    double angle = 0.5;

    double v2[] = {1.0,0.0,0.0,0.0 ,0.877582561890373  ,0.479425538604203,0.0 ,-0.479425538604203,0.877582561890373};
    Matrix m3 = R_x(angle);
    Matrix sol = Matrix(3, 3, v2, 9);

    printf("CREA\n");
    _assert(EqMatrix(m3,sol,10e-6) == true);

    return 0;
}

int RY()
{
    double angle = 0.5;

    double v2[] = {0.877582561890373,0.0,-0.479425538604203,0.0,1.0,0.0 ,0.479425538604203,0.0, 0.877582561890373};
    Matrix m3 = R_y(angle);
    Matrix sol = Matrix(3, 3, v2, 9);

    printf("CREA2\n");
    _assert(EqMatrix(m3,sol,10e-6) == true);

    return 0;
}

int RZ()
{
    double angle = 0.5;

    double v2[] = {0.877582561890373   ,0.479425538604203                   ,0,-0.479425538604203  , 0.877582561890373       ,            0,0                  , 0   ,1.000000000000000};
    Matrix m3 = R_z(angle);
    Matrix sol = Matrix(3, 3, v2, 9);

    printf("CREA2\n");
    _assert(EqMatrix(m3,sol,10e-6) == true);

    return 0;
}

int AccelPointMass(){

    double r_values[] = {1.0,2.0,3.0};
    double s_values[] = {1.0,1.0,1.0};

    double tolerance = 0.1;

    Matrix r = Matrix(1,3, r_values,3);
    Matrix s = Matrix(1,3, s_values, 3);

    Matrix result = AccelPointMass(r, s, Constants::GM_Earth);
    double num1 = -0.767106057663283e14;
    double num2 = -1.123624735995849e14;
    double num3 = -1.480143414328416e14;
    round(num3 * 100) / 100;
    round(num1 * 100) / 100;
    round(num2 * 100) / 100;
    double expected_values[] = {num1,num2,num3};
    Matrix expected(1,3, expected_values,3);

    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int AzELPa(){

    double r_values[] = {100,200,300};
    double tolerance = 10e-6;
    Matrix r = Matrix(1,3, r_values,3);

    Matrix result =  AzElPa(r);

    double res[] = {0.463647609000806,0.930274014115472,0.004000000000000,-0.002000000000000, 0.0,-0.000958314847500,  -0.001916629695000   ,0.001597191412500};
    Matrix expected(1,8, res,8);

    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int Cheb(){

    double t = 1.0;
    int N = 3;
    double Ta = 0.0;
    double Tb = 2.0;
    double tolerance = 10e-6;
    double Cx_values[] = {1.0, 2.0, 3.0};
    double Cy_values[] = {4.0, 5.0, 6.0};
    double Cz_values[] = {7.0, 8.0, 9.0};
    Matrix Cx(1, 3, Cx_values, 3);
    Matrix Cy(1, 3, Cy_values, 3);
    Matrix Cz(1, 3, Cz_values, 3);
    Matrix result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    double expected_values[] = {-2, -2.0, -2.0};
    Matrix expected(1, 3, expected_values, 3);
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int Legend(){


    double tolerance = 10e-6;
    Matrix result = Legendre(2,2,0.5);
    double expected_values[] = {1.000000000000000 , 0, 0,0.830389391308554 ,1.520017585030584  ,0,-0.347097518865836  , 1.629501555238869  , 1.491391294688047
            ,0         ,          0        ,           0,
            1.520017585030584,  -0.830389391308554     ,              0,
            2.822379484686224 ,  2.092581832544772 , -1.629501555238869};
    Matrix expected(1, 18, expected_values, 18);
    result.print();
    expected.print();
    printf("Legend");
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int EccA(){


    double tolerance = 10e-15;
    double expected = 1.934563210752024;
    double sol = EccAnom(1,1);
    _assert(fabs(expected-sol)<tolerance);

    return 0;
}

int MeanO(){


    double tolerance = 10e-15;
    double expected = 0.409413067075155;
    double sol = MeanObliquity(0.5);
    _assert(fabs(expected-sol)<tolerance);

    return 0;
}

int Geo(){


    double tolerance = 10e-4;
    double values[] = {1.0, 2.0, 3.0};
    Matrix r(1, 3, values, 3);
    Matrix result = Geodetic(r);
    double expected_values[] = {1.10714871779409, 1.570744136243924, -6.356748616533795e+06};
    Matrix expected(1, 3, expected_values, 3);
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}



int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);

    sol = m1 * m2;



    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));

    return 0;
}

int FracT()
{


    double x = 3.456;
    double sol = 0.456;
    double frac = Frac(x);
    double tol = 10e-3;
    _assert(fabs(frac-sol)<tol);

    return 0;
}

int Mjday1(){


    double tolerance = 10e-10;
    double expected =  5.267642372685205e+04;
    double sol = Mjday(2003,2,6,10,10,10);
    _assert(fabs(expected-sol)<tolerance);

    return 0;
}

int Mjday2(){


    double tolerance = 10e-10;
    double expected =  0.999999986586878;
    double sol = Mjday_TDB(1);
    _assert(fabs(expected-sol)<tolerance);

    return 0;
}

int Nut(){

    double tolerance = 10e-10;
    double d1 = 2.761314268655568e-05;
    double d2 =  3.961309941887382e-05;
    long double sol1;
    long double sol2;
    NutAngles(0.5,sol1,sol2);
    _assert((fabs(d1-sol1)<tolerance)&&(fabs(d2-sol2)<tolerance));

    return 0;
}

int Pos(){



    double tolerance = 10e-7;
    double v1[] = {1.0, 2.0, 3.0};
    double d1 = 1.866376490700467e+06 ;
    double d2 =  2.906709163731216e+06;
    double d3 =  5.343768700880938e+06;
    Position(1,1,1,v1);
    _assert((fabs(d1-v1[0])<tolerance)&&(fabs(d2-v1[1])<tolerance)&&(fabs(d3-v1[2])<tolerance));

    return 0;
}

int sign(){



    double d1 = -1;
    double sol = sign_(1,-3);
    _assert(d1==sol);

    return 0;
}

int timediffs(){


    double tolerance = 10e-7;
    double UT1_UTC = 1.0; // Provide your value
    double TAI_UTC = 1.0; // Provide your value
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    double sol1 = 0;
    double sol2 = 18;
    double sol3 = 19;
    double sol4 = 33.183999999999997;
    double sol5 = -18;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    _assert((fabs(UT1_TAI-sol1)<tolerance)&&(fabs(UTC_GPS-sol2)<tolerance)&&(fabs( UT1_GPS-sol3)<tolerance)&&(fabs(TT_UTC-sol4)<tolerance)&&(fabs( GPS_UTC-sol5)<tolerance));


    return 0;
}
int UNI(){


    double tolerance = 10e-4;
    double values[] = {1.0, 2.0, 3.0};
    Matrix r(1, 3, values, 3);
    Matrix result(1,3);
    result = unit(r);
    double expected_values[] = {0.267261241912424   ,0.534522483824849   ,0.801783725737273};
    Matrix expected(1, 3, expected_values, 3);
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}



int IER(){


    double tolerance = 10e-4;
    Matrix val(13,13);
    double Mjd_UTC = 37665.5;
    char interp = 'l';
    Matrix result(1,9);
    result = IERS(*globals::matrix,Mjd_UTC,interp);
    double expected_values[] = {-6.9328e-08,1.0353e-06,0.0323,0.0017,3.1086e-07,2.9954e-08,0.,0,2};
    Matrix expected(1, 9, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}
int all_tests()
{
    _verify(proMat_01);
    _verify(RX);
    _verify(RY);
    _verify(RZ);
    _verify(AccelPointMass);
    _verify(AzELPa);
    _verify(Cheb);
    _verify(FracT);
    _verify(Geo);
    _verify(EccA);
    _verify(MeanO);
    _verify(Mjday1);
    _verify(Mjday2);
    _verify(Nut);
    _verify(Pos);
    _verify(sign);
    _verify(timediffs);
    _verify(UNI);
    _verify(IER);
    _verify(Legend);
    return 0;
}


int main()
{
    globals::eop1962(13);


    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

   printf("Tests run: %d\n", tests_run);
   return result != 0;


}


