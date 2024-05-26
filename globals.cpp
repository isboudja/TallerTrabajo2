//
// Created by isboudja on 24/04/2024.
//

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "globals.h"
#include "SAT_Const.h"
#include "Mjday.h"

Matrix *globals::eopdata;
Matrix *globals::obs;
Matrix *globals::Rs;

void globals::eop1962(){
    globals::eopdata = new Matrix(13,21413);
    FILE *fid = fopen("../texts/eop19620101.txt", "r");

    if (fid == nullptr) {
        printf("error globals");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= 21413; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &globals::eopdata(1, i),
               &globals::eopdata(2, i), &globals::eopdata(3, i),
               &globals::eopdata(4, i), &globals::eopdata(5, i), &globals::eopdata(6, i),
               &globals::eopdata(7, i), &globals::eopdata(8, i), &globals::eopdata(9, i),
               &globals::eopdata(10, i), &globals::eopdata(11, i), &globals::eopdata(12, i),
               &globals::eopdata(13, i));
    }


    fclose(fid);

}

void globals::GEOS3(int nobs) {

    FILE *fid2 = fopen("../texts/GEOS3.txt", "r");
    globals::obs = new Matrix(nobs, 4);
    globals::Rs = new Matrix(nobs, 4);
    double y;
    double M;
    double d;
    double H;
    double min;
    double s;
    double az;
    double el;
    double Dist;
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    if (fid2 == nullptr) {
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    int ii;
    int i;
    int j;
    int l;
    for (n = 0; n < 46; n++) {
        fscanf(fid2, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &y, &M, &d, &H, &min, &s, &az, &el, &Dist);
        obs(n+1, 1) = Mjday(y, M, d, H, m, s);
        obs(n+1, 2) = Constants::Rad * az;
        obs(n+1, 3) = Constants::Rad * el;
        obs(n+1, 4) = 1e3 * Dist;

    }

    fclose(fid2);

}
/*
void globals::GGM(int c) {
//    = new Matrix(13,c);
FILE *fid = fopen("./texts/GGM03S.txt","r");
    globals::Cnm =  Cnm(181,181);
Matrix Snm(181,181);
Matrix temp(6,1);
if(fid==nullptr){
printf("error");
exit(EXIT_FAILURE);
}
int n;
int m;
for (n=0;n<=180;n++){
for (m=0;m<=n;m++) {
    fscanf(fid, "%d%d%lf%lf%lf%lf",&c,&c,&temp(1,3),&temp(1,4),&temp(1,5),&temp(1,6));
    Cnm(n+1,m+1) = temp(1,3);
    Snm(n+1,m+1) = temp(1,4);
}
}

fclose(fid);

}
 */

void globals::DE430() {
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
            fscanf(fid, "%lf",&PC(n+1,m+1));
        }
    }

    fclose(fid);

}
