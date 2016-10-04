#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include <R.h>
#include <time.h>
#include "functions.h"
//#ifndef M_PI
//#define M_PI 3.14159265358979323846
//#endif
//#define COMMENTS 0 

void EMstepFA(double *mu, double *alpha, double *sigmar, double *cpl, int *N, int *p, int *G, double *loglik, int *maxit, double *epsilon, int *label, double *pi, double *Lambdar, double *errr, int *q, int *v, double *x, int *counter){

      int i,j,g;
      int NN = *N;
      int GG = *G;
      int pp = *p;
      int qq = *q;
      double epsilonn= *epsilon;
      int maxitt = *maxit;
      int vv = *v;
 
      double *xx= malloc(sizeof(double)*NN*pp);
      double **sigma      = malloc(sizeof(double*)*GG);
      double **Lambda      = malloc(sizeof(double*)*GG);
      double **err      = malloc(sizeof(double*)*GG);
      for(g=0; g < GG; g++){
          sigma[g]        = malloc(sizeof(double)*pp*pp);
          Lambda[g]        = malloc(sizeof(double)*pp*qq);
          err[g]        = malloc(sizeof(double)*pp*pp);
       }
      for(i=0;i<NN;i++)
        for(j=0;j<pp;j++)
           xx[i*pp+j]= x[i+j*NN];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            sigma[g][i]=sigmar[g*pp*pp+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            err[g][i]=errr[g*pp*pp+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*qq; i++)
            Lambda[g][i]=Lambdar[g*pp*qq+i];
      i=2;
      while(getall(loglik,i)>epsilonn && i<(maxitt-1)){
           i=i+1;
           EMgrstepFA(xx,NN,pp,GG,vv,qq,label,mu,alpha,cpl,sigma,Lambda,err,pi);
           loglik[i]=llikFA(xx,mu,alpha,cpl,sigma,NN,pp,GG,pi);
      }
      *counter = i+1;
//      if(i<maxitt){
//         for(j=1;j<maxitt-i;j++)
//            loglik[i+j]=loglik[i];
//      }
     for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             sigmar[g*pp*pp+i]=sigma[g][i];
     for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             errr[g*pp*pp+i]=err[g][i];
     for(g=0; g<GG; g++)
         for(i=0; i<pp*qq; i++)
             Lambdar[g*pp*qq+i]=Lambda[g][i];
}

      void EMstep(double *mu, double *phi, double *alpha, double *cplr, double *gamr,
//                  double *sigmar, double *cpl0, double *wg, double *x, double *loglik,
                  double *cpl0, double *wg, double *x, double *loglik,
                  int *maxit, int *N, int *p, int *G, double *epsilon, int *label, int *v,
                  double *pi, int *counter){
      int i,j,g;
      int NN = *N;
      int GG = *G;
      int pp = *p;
      double epsilonn = *epsilon;
      int maxitt = *maxit;
      int vv = *v;
      double *xx= malloc(sizeof(double)*NN*pp);
      double **cpl      = malloc(sizeof(double*)*GG);
      double **gam      = malloc(sizeof(double*)*GG);
//      double **sigma      = malloc(sizeof(double*)*GG);
      for(g=0; g < GG; g++){
//          sigma[g]        = malloc(sizeof(double)*pp*pp);
          cpl[g]        = malloc(sizeof(double)*pp*2);
          gam[g]        = malloc(sizeof(double)*pp*pp);
      }

      for(i=0;i<NN;i++)
        for(j=0;j<pp;j++)
           xx[i*pp+j]= x[i+j*NN];
//      for(g=0; g<GG; g++)
//         for(i=0; i<pp*pp; i++)
//            sigma[g][i]=sigmar[g*pp*pp+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            gam[g][i]=gamr[g*pp*pp+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cpl[g][i]=cplr[g*pp*2+i];
      i=2;
      while(getall(loglik,i)>epsilonn && i<(maxitt-1)){
           i=i+1;
//           EMgrstep(xx,mu,alpha,cpl,cpl0,gam,sigma,wg,phi,NN,GG,pp,vv,pi,label,i);
           EMgrstep(xx,mu,alpha,cpl,cpl0,gam,wg,phi,NN,GG,pp,vv,pi,label,i+1);
//           loglik[i]=llik(xx,mu,alpha,phi,cpl0,cpl,gam,sigma,wg,NN,pp,GG,pi);
           loglik[i]=llik(xx,mu,alpha,phi,cpl0,cpl,gam,wg,NN,pp,GG,pi,label);
     }
      *counter=i+1; 
//      if(i<maxitt){
//         for(j=1;j<maxitt-i;j++)
//            loglik[i+j]=loglik[i];
//      }
     for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cplr[g*pp*2+i] = cpl[g][i];
     for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             gamr[g*pp*pp+i]=gam[g][i];
//     for(g=0; g<GG; g++)
//         for(i=0; i<pp*pp; i++)
//             sigmar[g*pp*pp+i]=sigma[g][i];

}

void EMstepSr(double *mu, double *phi, double *alpha, double *cplr, 
              double *gamr, double *x, double *loglik, int *maxit, int *N,
              int *p, int *G, double *eps, int *label, double *pi, int *v,
              int *counter){
      int i,g;
      int NN = *N;
      int GG = *G;
      int pp = *p;
      double epss = *eps;
      int maxitt = *maxit;
      int vv = *v;
      double **cpl      = malloc(sizeof(double*)*GG);
      double **gam      = malloc(sizeof(double*)*GG);
      for(g=0; g < GG; g++){
          cpl[g]        = malloc(sizeof(double)*pp*2);
          gam[g]        = malloc(sizeof(double)*pp*pp);
      }
      for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cpl[g][i]=cplr[g*pp*2+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            gam[g][i]=gamr[g*pp*pp+i];
      i=2;
      while(getall(loglik,i)>epss && i<(maxitt-1)){
           i=i+1;
           EMgrstepMSr(x,mu,phi,alpha,cpl,gam,NN,pp,GG,pi,vv,label,i+1);
           loglik[i]=llikMS(x,mu,alpha,cpl,gam,phi,NN,pp,GG,pi);
      }
      *counter = i+1;
//      if(i<maxitt){
//         for(j=1;j<maxitt-i;j++)
//            loglik[i+j]=loglik[i];
//      }
      for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cplr[g*pp*2+i] = cpl[g][i];

      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             gamr[g*pp*pp+i]=gam[g][i];

}
          
void EMstepMS(double *mu, double *phi, double *alpha, double *cplr,
              double *gamr, double *x, double *loglik, int *maxit, int *N,
              int *p, int *G, double *eps, int *label, double *pi, int *v,
              int *counter){
      int i,g;
      int NN = *N;
      int GG = *G;
      int pp = *p;
      double epss = *eps;
      int maxitt = *maxit;
      int vv = *v;
//      double *xx= malloc(sizeof(double)*NN*pp);
      double **cpl      = malloc(sizeof(double*)*GG);
      double **gam      = malloc(sizeof(double*)*GG);
      for(g=0; g < GG; g++){
          cpl[g]        = malloc(sizeof(double)*pp*2);
          gam[g]        = malloc(sizeof(double)*pp*pp);
      }
//      for(i=0;i<NN;i++)
//        for(j=0;j<pp;j++)
//           xx[i*pp+j]= x[i+j*NN];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cpl[g][i]=cplr[g*pp*2+i];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            gam[g][i]=gamr[g*pp*pp+i];
      i=2;
      while(getall(loglik,i)>epss && i<(maxitt-1)){
           i=i+1;
           EMgrstepMS(x,gam,mu,phi,alpha,cpl,NN,pp,GG,pi,vv,label,i+1);
           loglik[i]=llikMS(x,mu,alpha,cpl,gam,phi,NN,pp,GG,pi);
      }
      *counter = i+1;
//      if(i<maxitt){
//         for(j=1;j<maxitt-i;j++)
//            loglik[i+j]=loglik[i];
//      }
      for(g=0; g<GG; g++)
         for(i=0; i<pp*2; i++)
            cplr[g*pp*2+i] = cpl[g][i];

      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             gamr[g*pp*pp+i]=gam[g][i];
}

      void EMstepGH(double *mu, double *alpha, double *sigmar, double *cpl, int *N, int *p,
                   int *G, double *loglik, int *maxit, double *eps, int *label, double *pi,
                   int *v, double *x, int *counter){
      int i,j,g;
      int NN = *N;
      int GG = *G;
      int pp = *p;
      double epss = *eps;
      int maxitt = *maxit;
      int vv = *v;
      double *xx= malloc(sizeof(double)*NN*pp);
       double **sigma      = malloc(sizeof(double*)*GG);
      for(g=0; g < GG; g++)
          sigma[g]        = malloc(sizeof(double)*pp*pp);
      for(i=0;i<NN;i++)
        for(j=0;j<pp;j++)
           xx[i*pp+j]= x[i+j*NN];
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
            sigma[g][i]=sigmar[g*pp*pp+i];
      i=2;
      while(getall(loglik,i)>epss && i<(maxitt-1)){
           i=i+1;
           EMgrstepGH(xx,NN,pp,GG,vv,label,mu,alpha,cpl,sigma,pi);
           loglik[i]=llikGH(xx,mu,alpha,cpl,sigma,NN,pp,GG,pi);
      }
      *counter = i+1;
//      if(i<maxitt){
//         for(j=1;j<maxitt-i;j++)
//            loglik[i+j]=loglik[i];
///      }
      for(g=0; g<GG; g++)
         for(i=0; i<pp*pp; i++)
             sigmar[g*pp*pp+i]=sigma[g][i];

}








