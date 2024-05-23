//
// Created by isboudja on 22/05/2024.
//

#include <cstdio>
#include "AccelHarmonic.h"
#include "SAT_Const.h"
#include "Legendre.h"

Matrix AccelHarmonic(Matrix &r,Matrix &E,double n_max,double m_max) {

    FILE *fid = fopen("../texts/GGM03S.txt","r");
    Matrix Cnm(181,181);
    Matrix Snm(181,181);
    Matrix temp(1,6);
    Matrix a_bf(3,1);
    double c;
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

    Constants constants;

    double r_ref = 6378.1363e3;
    double gm = 398600.4415e9;
    Matrix r_bf =  r*E;
    r_bf.print();
    double d = r_bf.norm();
    d = pow(d,2.);
    double latgc = asin(r_bf(3,1) / d);
    double lon = atan2(r_bf(2,1), r_bf(1,1));

     Matrix Leg = Legendre(n_max, m_max, latgc);
     Matrix pnm(1,(n_max+1)*(m_max+1));
     Matrix dpnm(1,(n_max+1)*(m_max+1));
     for(int i=1;i<=(n_max+1)*(m_max+1);i++){
         pnm(1,i) = Leg(1,i);
     }
     int j = 1;
    for(int i=(n_max+1)*(m_max+1)+1;i<=((n_max+1)*(m_max+1))*2;i++){
        dpnm(1,j) = Leg(1,i);
        j++;
    }
    double dUdr = 0.;
    double dUdlatgc = 0.;
    double dUdlon = 0.;
    double q3 = 0.;
    double q2 = q3;
    double q1 = q2;
        double b1;
        j = 1;
    for(double n=0;n<=n_max;n++) {
        b1 = (-gm / pow(d,2)) * pow((r_ref / d),n ) * (n + 1);
        double b2 = (gm / d) * pow((r_ref / d),n);
        double b3 = (gm / d) * pow((r_ref / d),n);
        for(double m=0;m<=m_max;m++) {

            q1 = q1 + pnm(1, j) * (Cnm(n + 1, m + 1) * cos(m * lon) + Snm(n + 1, m + 1) * sin(m * lon));
            q2 = q2 + dpnm(1, j) * (Cnm(n + 1, m + 1) * cos(m * lon) + Snm(n + 1, m + 1) * sin(m * lon));
            q3 = q3 + m * pnm(1, j) * (Snm(n + 1, m + 1) * cos(m * lon) - Cnm(n + 1, m + 1) * sin(m * lon));
            j++;
        }
        dUdr = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon = dUdlon + q3 * b3;
        q3 = 0.;
        q2 = q3;
        q1 = q2;
    }

    double r2xy = pow(r_bf(1, 1),2) + pow(r_bf(2, 1),2);

    double ax = (1. / d * dUdr - r_bf(3,1) / (pow(d,2) * sqrt(r2xy)) * dUdlatgc) * r_bf(1,1) - (1. / r2xy * dUdlon) * r_bf(2,1);
    double ay = (1. / d * dUdr - r_bf(3,1) / (pow(d,2) * sqrt(r2xy)) * dUdlatgc) * r_bf(2,1) + (1. / r2xy * dUdlon) * r_bf(1,1);
    double az = 1. / d * dUdr * r_bf(3,1) + sqrt(r2xy) / pow(d,2) * dUdlatgc;

    a_bf(1,1) = ax;
    a_bf(2,1) = ay;
    a_bf(3,1) = az;


    }
    Matrix E2(1,3);
    for(int i=1;i<=3;i++){
        E2(1,i) = E(i,1);
    }
    a_bf.print();

    Matrix a = E2*a_bf;
    a.print();
    return a;
}



