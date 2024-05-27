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
#include "JPL_Eph_DE430.h"
#include "PoleMatrix.h"
#include "AccelHarmonic.h"
#include "Accel.h"
#include "DEInteg.h"
#include "G_AccelHarmonic.h"
#include "VarEqn.h"
#include "AuxParam.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "TimeUpdate.h"
#include "MeasUpdate.h"
#include "LTC.h"
#include "EqnEquinox.h"
#include "gast.h"
#include "GHAMatrix.h"
#include "gmst.h"

#define TOL_ 10e-14

using namespace std;
int tests_run = 0;
Constants c;

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

    Matrix result = AccelPointMass(r, s, c.GM_Earth);
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
    double sol1;
    double sol2;
    NutAngles(0.5,sol1,sol2);
    _assert((fabs(d1-sol1)<tolerance)&&(fabs(d2-sol2)<tolerance));

    return 0;
}

int Pos(){



    double tolerance = 10e-7;
    double v1[] = {1.0, 2.0, 3.0};
    Matrix r(1,3,v1,3);
    double d1 = 1.866376490700467e+06 ;
    double d2 =  2.906709163731216e+06;
    double d3 =  5.343768700880938e+06;
    Position(1,1,1,r);
    _assert((fabs(d1-r(1,1))<tolerance)&&(fabs(d2-r(1,2))<tolerance)&&(fabs(d3-r(1,3))<tolerance));

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
    result = IERS(*globals::eopdata,Mjd_UTC,interp);
    double expected_values[] = {-6.9328e-08,1.0353e-06,0.0323,0.0017,3.1086e-07,2.9954e-08,0.,0,2};
    Matrix expected(1, 9, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int PMatrix(){


    double tolerance = 10e-5;
    Matrix result = PrecMatrix(1,1);
    double expected_values[] = {1,0,0.,0.,1.,0,0.,0,1.};
    Matrix expected(3, 3, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;

}

int NMatrix(){


    double tolerance = 10e-4;
    Matrix val(13,13);
    double Mjd_UTC = 37665.5;
    char interp = 'l';
    Matrix result(3,3);
    result = NutMatrix(1.);
    double expected_values[] = {0.999999999623993     ,-2.51565103111378e-05   ,  -1.09162543714418e-05,
            2.51560791026998e-05       ,  0.999999998903467  ,    -3.9499840995938e-05,
            1.09172480376291e-05   ,    3.9499566370893e-05 ,        0.999999999160299};
    Matrix expected(3, 3, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int PoMatrix(){


    double tolerance = 10e-9;
    Matrix val(13,13);
    Matrix result = PoleMatrix(1.0,1.0);
    double expected_values[] = {0.54030230586814    ,     0.708073418273571       ,  0.454648713412841,
    0    ,      0.54030230586814    ,    -0.841470984807897,
    -0.841470984807897    ,     0.454648713412841   ,      0.291926581726429};
    Matrix expected(3, 3, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}
int JPL(){


    double tolerance = 10e2;
    double expected_values[] ={
            52891240655.1416, -82398828219.3078, -30762460775.264,
            34094359564.7165, -28358228564.6831, -9016395146.3542,
            -62008177373.3339, 122716094121.445, 53207795594.9921,
            -160126907125.025, -22381470230.0169, -1176302999.30631,
            583165485249.939, -624231961619.039, -280912964474.469,
            -1289046343232.77, 191634032223.748, 134639477339.296,
            -131859240183.77, 2466292748344.82, 1083466135600.11,
            -4288263038254.25, -1332516416256.53, -440122288690.766,
            -3909996317326.24, 2903217225254.62, 2087565919619.75,
            -115573196.913284, -310886551.21669, -165562926.423869,
            62129751075.8268, -122378505237.373, -53074092053.2144
    };
    Matrix expected(1, 33, expected_values, 33);
    expected.print();
    Matrix result(1,33);
    result = JPL_Eph_DE430(33264.0);
    result.print();
    //-2.1528481871834e+23
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int GHMC(){


    double tolerance = 10e-6;
    double E2[] = {

                          -0.984320311914274,         0.176389708356208   ,  -0.000440838970817921,
    -0.176389673454202,         -0.98432040991552   ,   -0.00011714291893872,
    -0.000454589601769828 ,    -3.75467123998344e-05    ,     0.999999895969264};
    double expected_values2[] ={
            1,1,1};
    Matrix r(3,1,expected_values2,3);
    Matrix E(3,3,E2,9);
    Matrix result(3,3);
    result = G_AccelHarmonic(r,E,2,2);
    result.print();
    _assert(EqMatrix(result, E, tolerance));

    return 0;
}



int VEQ(){


    double tolerance = 10e-6;
    double expected_values[] ={


            1,
            1,
            1,
            2.6346498390509e+138,
            -3.23987314882573e+139,
            -6.97967557748611e+139 ,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140,
            1,
            1,
            1,
            8.88176580222003e+140,
            -1.78113267704466e+140,
            6.04647582287588e+140};

    double resul[] ={


            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,};
    Matrix expected(42, 1, expected_values, 42);
    Matrix Y(42, 1, resul, 42);
    Matrix result(42,1);
    result = VarEqn(1,Y);
    expected.print();
    result.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}


int HMC(){


    double tolerance = 10e-6;
    double r2  [] ={7101597.83995656,
            1295244.79268407,
            12755.6074812101};
    double expected_values2[] ={   -7.53454836070451,
                                   -1.37427214614486,
                                   -0.0135587742125453};
    double E2[] ={    -0.984333171840543    ,     0.176317930080741 ,    -0.000440847341518802,
    -0.17631789518587     ,   -0.984333269844254     ,-0.000117110813026442,
    -0.000454589441322244    , -3.75467826882455e-05   ,      0.999999895969334};
    Matrix r(3,1,r2,3);
    Matrix E(3,3,E2,9);
    Matrix expected(1,3,expected_values2,3);
    Matrix result = AccelHarmonic(r,E,2,2);
    expected.print();
    result.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int ASCE(){


    double tolerance = 10e-6;
    double Y2[] ={1,1,1,1,1,1};
    double expected_values[] ={  1,
            1,
            1,
            2.59173969179228e+138,
            -3.23976995737588e+139,
            -6.9808470022172e+139};
    Matrix expected(1, 6, expected_values, 6);
    Matrix Y(6,1,Y2,6);
    Matrix res=Accel(1, Y);
    expected.print();
    res.print();
    _assert(EqMatrix(res, expected, tolerance));

    return 0;
}

int DEI(){


    double tolerance = 10e-6;
    double Y2[] ={6221397.62857869,2867713.77965738,3006155.98509949,4645.04725161807,-2752.21591588205,-7507.99940987033};
    double expected_values[] ={  5542555.89427452,
            3213514.83814162,
            3990892.92789074,
            5394.0689404439,
            -2365.21290574022,
            -7061.84481373472};
    Matrix expected(6, 1, expected_values, 6);
    Matrix Y(6,1,Y2,6);

    DEInteg(Accel,0,-134.999991953373,1e-13,1e-6,6,Y);
    expected.print();
    Y.print();
    _assert(EqMatrix(Y, expected, tolerance));

    return 0;
}

int Tup(){


    double tolerance = 10e-9;
    Matrix P(1,1);
    P(1,1)  = 1.;
    Matrix Phi(1,1);
    Phi(1,1) = 1.;
    Matrix result = TimeUpdate(P,Phi,1.0);
    double expected_values[] = {2};
    Matrix expected(1, 1, expected_values, 1);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int Mup(){


    double tolerance = 10e-9;
    Matrix x(1,1);
    x(1,1)  = 1.;
    Matrix P(1,1);
    P(1,1) = 1.;
    Matrix g(1,1);
    g(1,1)  = 1.;
    Matrix G(1,1);
    G(1,1) = 1.;
    Matrix result = MeasUpdate(x,1.,g,1.,G,P,1.);
    double expected_values[] = {0.5};
    Matrix expected(1, 1, expected_values, 1);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int LC(){


    double tolerance = 10e-9;
    Matrix result = LTC(1.,1.);
    double expected_values[] = {-0.841470984807897  ,        0.54030230586814    ,                     0,
    -0.454648713412841     ,   -0.708073418273571 ,         0.54030230586814,
    0.291926581726429      ,   0.454648713412841      ,   0.841470984807897};
    Matrix expected(1, 1, expected_values, 1);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;
}

int Nox(){


    double tolerance = 10e-9;
    double result = EqnEquinox(1.);
    double expected_values[] = {2.51565103142908e-05};
    Matrix res(1,1);
    res(1,1) = result;
    Matrix expected(1, 1, expected_values, 1);
    res.print();
    expected.print();
    _assert(EqMatrix(res, expected, tolerance));

    return 0;
}

int GAS(){

    double tolerance = 10e-9;
    double result = gast(1.);
    double expected_values[] = {0.990436095961894};
    Matrix res(1,1);
    res(1,1) = result;
        Matrix expected(1, 1, expected_values, 1);
        res.print();
        expected.print();
        _assert(EqMatrix(res, expected, tolerance));

        return 0;

        }

int GHAM(){

    double tolerance = 10e-9;
    Matrix result = GHAMatrix(1.);
    double expected_values[] = {0.548325220865005        , 0.836265180527889        ,                 0,
    -0.836265180527889  ,       0.548325220865005   ,                      0,
    0            ,             0    ,                     1};
    Matrix expected(3, 3, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;

}

int GBS(){

    double tolerance = 10e-9;
    Matrix result = GHAMatrix(1.);
    double expected_values[] = {0.548325220865005        , 0.836265180527889        ,                 0,
                                -0.836265180527889  ,       0.548325220865005   ,                      0,
                                0            ,             0    ,                     1};
    Matrix expected(3, 3, expected_values, 9);
    result.print();
    expected.print();
    _assert(EqMatrix(result, expected, tolerance));

    return 0;

}

int gms(){

    double tolerance = 10e-9;
    double result = gmst(1.);
    double expected_values[] = { 0.99041093945158};
    Matrix res(1,1);
    res(1,1) = result;
    Matrix expected(1, 1, expected_values, 1);
    res.print();
    expected.print();
    _assert(EqMatrix(res, expected, tolerance));

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
    _verify(PMatrix);
    _verify(PoMatrix);
    _verify(NMatrix);
    _verify(Tup);
    _verify(Mup);
    _verify(Nox);
    _verify(GAS);
    _verify(GHAM);
    _verify(GBS);
    _verify(gms);
    _verify(HMC);
//    _verify(VEQ);
 //   _verify(GHMC);
   // _verify(ASCE);
   // _verify(DEI);
    return 0;
}


int main()
{
    AuxParam AuxParam;
    globals::eop1962();
    globals::GGM();
    globals::DE430();
    globals::GEOS3(46);

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

   printf("Tests run: %d\n", tests_run);
   return result != 0;

}


