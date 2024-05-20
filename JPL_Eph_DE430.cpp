//
// Created by isboudja on 15/05/2024.
//

#include <cstdio>
#include <cstdlib>
#include "JPL_Eph_DE430.h"
#include "Cheb3D.h"


Matrix JPL_Eph_DE430(double Mjd_TDB) {
    // Definición del vector de resultado
   Matrix result(11, 1.0); // r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun

    FILE *fid = fopen("./texts/DE430Coeff.txt","r");
    Matrix PC(2285,1020);
    if(fid==nullptr){
        printf("error2");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n=0;n<=2285;n++){
        for (m=0;m<=1020;m++) {
            fscanf(fid, "%lf",&PC(n,m));
        }
    }

    fclose(fid);

    // Obtener el día juliano
    double JD = Mjd_TDB + 2400000.5;
    int i;
    for (i = 0; i < 13; i++) {
        if(PC(i+1,1)<=JD && JD<=PC(i+1,2)){
            break;

        }
    }

    int j;
    Matrix PCtemp(1,1020);
    for (j = 0; i < 1020; i++) {
         PCtemp(1,i+1) = PC(i,j+1);
    }
      ;


    double t1 = PCtemp(1,1)-2400000.5;


    double  dt = Mjd_TDB - t1;
    Matrix temp(1,4);
    n=231;
    j=1;
    while(n<=270){
        temp(1,j) = n;
        j++;
        n += 13;
    }
    Matrix Cx_Earth=    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Earth = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Earth = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    for(j=1;j<temp.col;j++){
        temp(1,j) += 39;
    }

    Matrix Cx= PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz= PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);


        Matrix x = Matrix::concat(Cx_Earth,Cx);
        Matrix y  = Matrix::concat(Cy_Earth,Cz);
        Matrix z = Matrix::concat(Cz_Earth,Cy);

    double Mjd0;
    if (0<=dt && dt<=16) {
        j = 0;
         Mjd0 = t1;
    }else {
        if (16 < dt && dt <= 32) {
            j = 1;
            Mjd0 = t1 + 16 * j;
        }
    }
    Matrix r_Earth(1,3);
    Matrix res1 = x.sub((13*j+1),(13*j+13));
    Matrix res2= y.sub((13*j+1),(13*j+13));
    Matrix res3= z.sub((13*j+1),(13*j+13));

    r_Earth =Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, res1,res2,res3)*1e3;

    n=441;
    j=1;
    while(n<=480){
        temp(1,j) = n;
        j++;
        n += 13;
    }
    Matrix Cx_Moon=    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Moon = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Moon = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);



    for(i=1;i<=7;i++){

        for(j=1;j<temp.col;j++){
            temp(1,j) += 39;
        }
         Cx= PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
         Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
         Cz= PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
         x = Matrix::concat(Cx_Moon,Cx);
         y  = Matrix::concat(Cy_Moon,Cz);
         z = Matrix::concat(Cz_Moon,Cy);
    }

    if (0<=dt && dt<=4) {
        j = 0;
        Mjd0 = t1;
    }else {
        if (4 < dt && dt <= 8) {

            j = 1;
            Mjd0 = t1 + 4 * j;
        } else {
            if (8 < dt && dt <= 12) {
                j = 2;
                Mjd0 = t1 + 4 * j;
            } else {
                if (12 < dt && dt <= 16) {


                    j = 3;
                    Mjd0 = t1 + 4 * j;
                } else {
                    if (16 < dt && dt <= 20) {
                        j = 4;
                        Mjd0 = t1 + 4 * j;
                    } else {
                        if (20 < dt && dt <= 24) {
                            j = 5;
                            Mjd0 = t1 + 4 * j;
                        } else {
                            if (24 < dt && dt <= 28) {
                                j = 6;
                                Mjd0 = t1 + 4 * j;
                            } else {
                                if (28 < dt && dt <= 32) {
                                    j = 7;
                                    Mjd0 = t1 + 4 * j;
                                }
                            }
                        }
                    }
                }
            }
        }
    }





    Matrix r_Moon(1,3);
    res1 = x.sub((13*j+1),(13*j+13));
     res2= y.sub((13*j+1),(13*j+13));
     res3= z.sub((13*j+1),(13*j+13));
    r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, res1,res2,res3)*1e3;


    n=753;
    j=1;
    while(n<=786){
        temp(1,j) = n;
        j++;
        n += 11;
    }

    Matrix Cx_Sun=    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Sun = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Sun = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    for(j=1;j<temp.col;j++){
        temp(1,j) += 39;
    }
    Cx=PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Cz= PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    x = Matrix::concat(Cx_Sun,Cx);
    y  = Matrix::concat(Cy_Sun,Cz);
    z = Matrix::concat(Cz_Sun,Cy);
    if (0<=dt && dt<=16) {
        j = 0;
        Mjd0 = t1;
    }else {
        if (16 < dt && dt <= 32) {
            j = 1;
            Mjd0 = t1 + 16 * j;
        }
    }
    Matrix r_Sun(1,3);
    r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, res1,res2,res3)*1e3;

    n=3;
    j=1;
    while(n<=45){
        temp(1,j) = n;
        j++;
        n += 14;
    }
    Matrix Cx_Mercury =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Mercury  = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Mercury  = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    for(i=1;i<=3;i++){

        for(j=1;j<temp.col;j++){
            temp(1,j) += 42;
        }
        Cx= PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
        Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
        Cz= PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
        x = Matrix::concat(Cx_Mercury,Cx);
        y  = Matrix::concat(Cy_Mercury,Cz);
        z = Matrix::concat(Cz_Mercury,Cy);
    }

    if (0<=dt && dt<=8){
        j=0;
    Mjd0 = t1;
    }
    else{if(8<dt && dt<=16){
    j=1;
    Mjd0 = t1+8*j;
    }
    else{if (16<dt && dt<=24){
    j=2;
    Mjd0 = t1+8*j;
                }
    else{if(24<dt && dt<=32){
    j=3;
    Mjd0 = t1+8*j;
                    }
                }
        }
    }
    Matrix r_Mercury(1,3);
    res1 = x.sub((14*j+1),(14*j+14));
    res2= y.sub((14*j+1),(14*j+14));
    res3= z.sub((14*j+1),(14*j+14));
    r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, res1,res2,res3)*1e3;

    n=171;
    j=1;
    while(n<=201){
        temp(1,j) = n;
        j++;
        n += 10;
    }
    Matrix Cx_Venus =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Venus = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Venus = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    for(j=1;j<temp.col;j++){
        temp(1,j) += 30;
    }
    Cx= PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Cz= PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    x = Matrix::concat(Cx_Venus,Cx);
    y  = Matrix::concat(Cy_Venus,Cz);
    z = Matrix::concat(Cz_Venus,Cy);
    if (0<=dt && dt<=16){
        j=0;
    Mjd0 = t1;
        }
    else{if(16<dt && dt<=32){
    j=1;
    Mjd0 = t1+16*j;
    }
            }
    Matrix r_Venus(1,3);
    res1 = x.sub((10*j+1),(10*j+10));
    res2= y.sub((10*j+1),(10*j+10));
    res3= z.sub((10*j+1),(10*j+10));
    r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, res1,res2,res3)*1e3;


    n=309;
    j=1;
    while(n<=342){
        temp(1,j) = n;
        j++;
        n += 11;
    }
    Matrix Cx_Mars =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Mars =PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Mars = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Mars(1,3);
    res1 = Cx_Mars.sub((11*j+1),(11*j+11));
    res2= Cy_Mars.sub((11*j+1),(11*j+11));
    res3= Cz_Mars.sub((11*j+1),(11*j+11));
    r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, res1,res2,res3)*1e3;

    n=342;
    j=1;
    while(n<=366){
        temp(1,j) = n;
        j++;
        n += 8;
    }
    Matrix Cx_Jupiter =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Jupiter = PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Jupiter = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Jupiter(1,3);
    res1 = Cx_Jupiter.sub((8*j+1),(8*j+8));
    res2= Cy_Jupiter.sub((8*j+1),(8*j+8));
    res3= Cz_Jupiter.sub((8*j+1),(8*j+8));
    r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32,res1,res2,res3)*1e3;

    n=366;
    j=1;
    while(n<=387){
        temp(1,j) = n;
        j++;
        n += 7;
    }

    Matrix Cx_Saturn =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Saturn= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Saturn = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Saturn(1,3);
    res1 = Cx_Saturn.sub((7*j+1),(7*j+7));
    res2= Cy_Saturn.sub((7*j+1),(7*j+7));
    res3= Cz_Saturn.sub((7*j+1),(7*j+7));
    r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32,res1,res2,res3)*1e3;


    n=387;
    j=1;
    while(n<=405){
        temp(1,j) = n;
        j++;
        n += 6;
    }
    Matrix Cx_Uranus =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Uranus= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Uranus = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Uranus(1,3);
    res1 = Cx_Uranus.sub((6*j+1),(6*j+6));
    res2= Cy_Uranus.sub((6*j+1),(6*j+6));
    res3= Cz_Uranus.sub((6*j+1),(6*j+6));
    r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, res1,res2,res3)*1e3;

    n=405;
    j=1;
    while(n<=423){
        temp(1,j) = n;
        j++;
        n += 6;
    }
    Matrix Cx_Neptune =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Neptune= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Neptune = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Neptune(1,3);
    res1 = Cx_Neptune.sub((6*j+1),(6*j+6));
    res2= Cy_Neptune.sub((6*j+1),(6*j+6));
    res3= Cz_Neptune.sub((6*j+1),(6*j+6));
    r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,res1,res2,res3)*1e3;

    n=423;
    j=1;
    while(n<=441){
        temp(1,j) = n;
        j++;
        n += 6;
    }

    Matrix Cx_Pluto =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Pluto= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Pluto = PCtemp.sub((int)temp(1,4),(int)temp(1,3)-1);
    j=0;
    Mjd0 = t1;
    Matrix r_Pluto(1,3);
    res1 = Cx_Pluto.sub((6*j+1),(6*j+6));
    res2= Cy_Pluto.sub((6*j+1),(6*j+6));
    res3= Cz_Pluto.sub((6*j+1),(6*j+6));
    r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,res1,res2,res3)*1e3;

    n=819;
    j=1;
    while(n<=839){
        temp(1,j) = n;
        j++;
        n += 10;
    }

    Matrix Cx_Nutations =    PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Nutations= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
    for(i=1;i<=3;i++){
        for(j=1;j<temp.col;j++){
            temp(1,j) += 20;
        }
        Cx= PCtemp.sub((int)temp(1,2),(int)temp(1,1)-1);
        Cy= PCtemp.sub((int)temp(1,3),(int)temp(1,2)-1);
        x = Matrix::concat(Cx_Nutations,Cx);
        y  = Matrix::concat(Cy_Nutations,Cz);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }
    else{if(8<dt && dt<=16){
            j=1;
            Mjd0 = t1+8*j;
        }
        else{if (16<dt && dt<=24){
                j=2;
                Mjd0 = t1+8*j;
            }
            else{if(24<dt && dt<=32){
                    j=3;
                    Mjd0 = t1+8*j;
                }
            }
        }
    }

    Matrix Nutations(1,3);
    res1 = Cx_Pluto.sub((10*j+1),(10*j+10));
    res2= Cy_Pluto.sub((10*j+1),(10*j+10));
    Matrix res32(1,10);
    Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, res1,res2,res32);



    double EMRAT = 81.30056907419062;
    double EMRAT1 = 1. / (1. + EMRAT);

    r_Earth = r_Earth - (r_Moon*EMRAT1);
    r_Mercury =  r_Mercury-r_Earth ;
    r_Venus =   r_Venus-r_Earth;
    r_Mars =   r_Mars-r_Earth;
    r_Jupiter =  r_Jupiter -r_Earth;
    r_Saturn =   r_Saturn-r_Earth;
    r_Uranus =   r_Uranus-r_Earth;
    r_Neptune =  r_Neptune-r_Earth;
    r_Pluto =   r_Pluto-r_Earth;
    r_Sun =   r_Sun-r_Earth;

    Matrix result2(1,27);
    result2 = Matrix::concat(r_Earth,r_Mercury);
    result2 = Matrix::concat(result2,r_Venus);
    result2 = Matrix::concat(result2,r_Mars);
    result2 = Matrix::concat(result2,r_Jupiter);
    result2 = Matrix::concat(result2,r_Saturn);
    result2 = Matrix::concat(result2,r_Uranus);
    result2 = Matrix::concat(result2,r_Neptune);
    result2 = Matrix::concat(result2,r_Pluto);
    result2 = Matrix::concat(result2,r_Sun);


    return result2;
}