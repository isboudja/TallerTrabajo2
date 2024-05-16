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
        printf("error");
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
    Matrix Cx_Earth=    temp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Earth = temp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Earth = temp.sub((int)temp(1,4),(int)temp(1,3)-1);
    temp(1,j+1) += 39;
    Matrix Cx= temp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy= temp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz= temp.sub((int)temp(1,4),(int)temp(1,3)-1);


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
    Matrix Cx_Moon=    temp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Moon = temp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Moon = temp.sub((int)temp(1,4),(int)temp(1,3)-1);



    for(i=1;i<=7;i++){

        temp(1,j+1) += 39;
         Cx= temp.sub((int)temp(1,2),(int)temp(1,1)-1);
         Cy= temp.sub((int)temp(1,3),(int)temp(1,2)-1);
         Cz= temp.sub((int)temp(1,4),(int)temp(1,3)-1);
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

    Matrix Cx_Sun=    temp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Matrix Cy_Sun = temp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Matrix Cz_Sun = temp.sub((int)temp(1,4),(int)temp(1,3)-1);
    temp(1,j+1) += 39;
    Cx= temp.sub((int)temp(1,2),(int)temp(1,1)-1);
    Cy= temp.sub((int)temp(1,3),(int)temp(1,2)-1);
    Cz= temp.sub((int)temp(1,4),(int)temp(1,3)-1);
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

    // Calcular r_Mercury usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Venus usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Mars usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Jupiter usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Saturn usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Uranus usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Neptune usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Pluto usando Cheb3D y los coeficientes correspondientes

    // Aplicar correcciones EMRAT

    return result;
}