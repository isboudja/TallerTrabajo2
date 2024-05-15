//
// Created by isboudja on 15/05/2024.
//

#include <cstdio>
#include <cstdlib>
#include "JPL_Eph_DE430.h"


Matrix JPL_Eph_DE430(double Mjd_TDB) {
    // Definición del vector de resultado
   Matrix result(11, 1.0); // r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun

    FILE *fid = fopen("./texts/DE430","r");
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
    Matrix Cx_Earth(1,4);
    Matrix Cy_Earth(1,4);
    Matrix Cz_Earth(1,4);
    Matrix Cx(1,2);
    Matrix Cy(1,2);
    Matrix Cz(1,2);

    for (j = 0; j < 2; j++) {
        Cx_Earth(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 1; j < 3; j++) {
        Cy_Earth(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 2; j < 4; j++) {
        Cz_Earth(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }

    for (j=0;j<temp.col;j++){
        temp(1,j+1) += 39;
    }

    for (j = 0; j < 2; j++) {
        Cx(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 1; j < 3; j++) {
        Cy(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 2; j < 4; j++) {
        Cz(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }

    for (j = 2; j < 4; j++) {
        Cx_Earth(1,j+1)  = Cx(1,j-1);
        Cy_Earth(1,j+1)  = Cy(1,j-1);
        Cz_Earth(1,j+1)  = Cz(1,j-1);
    }
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
    r_Earth = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth(13*j+1:13*j+13),Cy_Earth(13*j+1:13*j+13), Cz_Earth(13*j+1:13*j+13));

    n=441;
    j=1;
    while(n<=480){
        temp(1,j) = n;
        j++;
        n += 13;
    }
    Matrix Cx_Moon(1,4);
    Matrix Cy_Moon(1,4);
    Matrix Cz_Moon(1,4);


    for (j = 0; j < 2; j++) {
        Cx_Moon(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 1; j < 3; j++) {
        Cy_Moon(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 2; j < 4; j++) {
        Cz_Moon(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for(i=1;i<=7;i++){
    for (j=0;j<temp.col;j++){
        temp(1,j+1) += 39;
    }

    for (j = 0; j < 2; j++) {
        Cx(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 1; j < 3; j++) {
        Cy(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }
    for (j = 2; j < 4; j++) {
        Cz(1,j+1)  = PCtemp(1,temp(1,j+1)-1);
    }

    for (j = 2; j < 4; j++) {
        Cx_Moon(1,j+1)  = Cx(1,j-1);
        Cy_Moon(1,j+1)  = Cy(1,j-1);
        Cz_Moon(1,j+1)  = Cz(1,j-1);
    }
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
    r_Moon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon(13*j+1:13*j+13),Cy_Moon(13*j+1:13*j+13), Cz_Moon(13*j+1:13*j+13));

    // Calcular r_Moon usando Cheb3D y los coeficientes correspondientes

    // Calcular r_Sun usando Cheb3D y los coeficientes correspondientes

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