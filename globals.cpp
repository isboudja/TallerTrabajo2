//
// Created by isboudja on 24/04/2024.
//

#include <cstdio>
#include <cstdlib>
#include "globals.h"

Matrix *globals::matrix;
Matrix *globals::matrix2;

void globals::eop1962(int c) {
    globals::matrix  = new Matrix(13,c);
    FILE *fid = fopen("D:/UNI/2do Semestre/TTT/TrabajoTaller/texts/eop19620101.txt","r");

            if(fid==nullptr){
                printf("error");
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



void globals::egm(int c) {
//    globals::matrix2  = new Matrix(13,c);
    FILE *fid = fopen("egm.txt","r");
    Matrix obs(10,4),Rs(3,1);
    double Cnm[352][362];
    if(fid==nullptr){
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n=0;n<=180;n++){
        for (m=0;m<=n;m++) {
            fscanf(fid, "%d%d%lf%lf%lf%lf",&c,&c,&Cnm[n+1][m+1]);
        }
    }

    fclose(fid);

}