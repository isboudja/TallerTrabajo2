//
// Created by isboudja on 24/04/2024.
//

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "globals.h"
#include "SAT_Const.h"
#include "Mjday.h"

Matrix *globals::matrix;
Matrix *globals::matrix2;

void globals::eop1962(int c) {
    globals::matrix  = new Matrix(13,c);
    FILE *fid = fopen("../texts/eop19620101.txt","r");

            if(fid==nullptr){
                printf("error globals");
                exit(EXIT_FAILURE);
            }
            for(int i=1;i<=c;i++){
                fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(*globals::matrix)(1,i),&(*globals::matrix)(2,i),&(*globals::matrix)(3,i),
                                                                                                                                   &(*globals::matrix)(4,i),&(*globals::matrix)(5,i),&(*globals::matrix)(6,i)
                                                                                                                                                                                    ,&(*globals::matrix)(7,i),&(*globals::matrix)(8,i),&(*globals::matrix)(9,i)
                                                                                                                                                                                                                                        ,&(*globals::matrix)(10,i),&(*globals::matrix)(11,i),&(*globals::matrix)(12,i),&(*globals::matrix)(13,i));
            }


    fclose(fid);

}



void globals::GEOS3(int c) {
    FILE *fid = fopen("../texts/GEOS3.txt","r");
    Matrix obs(c,4),Rs(1,3);
    double Y;
    double M;
    double d;
    double H;
    double min;
    double s;
    double az;
    double el;
    double Dist;
    double Cnm[352][362];
    if(fid==nullptr){
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n=0;n<=46;n++){
            fscanf(fid, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",&Y,&M,&d,&H,&min,&s,&az,&el,%Dist);
            obs(n,1) = Mjday(Y,M,d,H,m,s);
            obs(n,2) = Constants::Rad*az;
            obs(n,3) = Constants::Rad*el;
            obs(n,4) = 1e3*Dist;

    }

    fclose(fid);

}

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