#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include <R.h>
#include <time.h>
#include "functions.h"
#include "zbsubs.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#endif
#define COMMENTS 0 



      void BesselK(double zr, int nz, double nu, int expon, int nSeq, double *r, double *r1){
      double *zi= malloc(sizeof(double)*nz);
      double *cyr= malloc(sizeof(double)*nSeq);
      double *cyi= malloc(sizeof(double)*nSeq);
      int i,j;
      int ierr =1;
      int n = nSeq;
      int kode = expon+(double)1.0;
      if(nu <0.0) nu = -nu;
      double fnu = nu;
      for(i=0;i<nz;i++){
          zi[i]=0.0;
          r[i]=0.0;
       }
       for(i=0;i<nSeq;i++)
          cyr[i]=0.0;
       for(i=0;i<nSeq;i++)
          cyi[i]=0.0;
       for(i=0;i<nz;i++){
          zbesk_(&zr,&zi[i],&fnu,&kode,&n,cyr,cyi,&nz,&ierr);
          if(nSeq >1){
            for(j=0;j<nSeq;j++){
               r1[i*nSeq+j]=cyr[j];
            }
          }
          else
          r[i]=cyr[0];
        }

      free(zi); free(cyr); free(cyi);
}

      void BesselK1(double *zr, int nz, double nu, int expon, int nSeq, double *r, double *r1){
      double *zi= malloc(sizeof(double)*nz);
      double *cyr= malloc(sizeof(double)*nSeq);
      double *cyi= malloc(sizeof(double)*nSeq);
      int i,j;
      int ierr =1;
      int n = nSeq;
      int kode = expon+(double)1.0;
      if(nu <0.0) nu = -nu;
      double fnu = nu;
      for(i=0;i<nz;i++){
          zi[i]=0.0;
          r[i]=0.0;
       }
       for(i=0;i<nSeq;i++)
          cyr[i]=0.0;
       for(i=0;i<nSeq;i++)
          cyi[i]=0.0;
       for(i=0;i<nz;i++){
          zbesk_(&zr[i],&zi[i],&fnu,&kode,&n,cyr,cyi,&nz,&ierr);
          if(nSeq >1){
            for(j=0;j<nSeq;j++){
               r1[i*nSeq+j]=cyr[j];
            }
          }
          else
          r[i]=cyr[0];
        }
       free(zi); free(cyr); free(cyi);
}



void cov(double *x, double *cov, int N, int p){
      int j,k,i;
      double *wt = (double*)malloc(sizeof(double)*p);
     for(i=0; i<p; i++){
        wt[i]=0;
        for(j=0; j<N; j++){
            wt[i] += x[i+j*p];
        }
        wt[i] /=N;
     }
     for(j=0; j<p; j++){
        for(k=0; k<p; k++){
           for(i=0; i<N; i++){
           cov[j*p+k] += (x[j+i*p]-wt[j])*(x[k+i*p]-wt[k]);
           }
           cov[j*p+k] /=((double)N-(double)1.0);
        }
     }
     free(wt);
}
double getall(double loglik[], int n) {
//      if(n < 2) exit(0);
      double lm1 = loglik[n];
      double lm  = loglik[n-1];
      double lm_1= loglik[n-2];
      double am  = (lm1 - lm)/(lm - lm_1);
      double lm1Inf = lm + (lm1 - lm)/(1-am);
      double val =  lm1Inf - lm;
      if(val == NAN) val = 0;
      if(val < 0) val =1;
      return(val);
}
/*void rgpar(double *x, int N, int p, int G, double *w, int q){
      int i,g;
      double *mu = malloc(sizeof(double)*G*p);
      double *mu0 = malloc(sizeof(double)*p);
      double *pi = malloc(sizeof(double)*G);
      double *wt = malloc(sizeof(double)*N);
      double *cpl = malloc(sizeof(double)*G*2);
      double *alpha = malloc(sizeof(double)*G*p);
      double **err    = malloc(sizeof(double*)*G);
      double **sigma    = malloc(sizeof(double*)*G);
      double **Lambda    = malloc(sizeof(double*)*G);
      for(g=0; g<G; g++){
         err[g]      = malloc(sizeof(double)*p*p);
         sigma[g]      = malloc(sizeof(double)*p*p);
         Lambda[g]      = malloc(sizeof(double)*p*q);
      }

      if(!w){
        for(i=0; i<N; i++)
           for(g=0; g<G; g++)
              w[i*G+g]= (double)1.0/(double)G;     
      }
      for(g=0; g<G; g++){
        for(i=0; i<N; i++){
           wt[i]=w[g+i*G];
        }
      iparM(x,N,p,wt,q,G,mu0,alpha,Lambda[g],sigma[g],err[g],cpl);
      for(i=0;i<p;i++){
         mu[g*p+i]=mu[i];
      }
}
      for(g=0; g<G; g++)
         pi[g]=(double)1.0/(double)G;
}*/


/* double randn (double mu, double sigma)

{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;
  return (mu + sigma * (double) X1);
}*/

void weightedmean(double *x, int N, int p, double *wt, double *wmean){
      int i,j;
      double sumwt;
      for(j =0; j<p; j++){
      sumwt =0.0;
          wmean[j] = 0.0;
          for(i=0; i<N; i++){
//          wmean[j] += x[j+i*p]*wt[j];
          wmean[j] += x[j+i*p]*wt[i];
//          sumwt = sumwt + wt[j];
          sumwt = sumwt + wt[i];
          }
      wmean[j]/= sumwt;
      }
}
//void Covariance(int N, int p, int G, double *x, double *z, double *mu, int g, double *sampcov) {
//  int j,k,i;
//        double sum;
//        double *Wt = (double*)malloc(sizeof(double)*N);
//
//  sum=0;
//        for(i=0; i < N; i++) {
//                Wt[i] = z[i + N*g];
//                sum  += Wt[i];
//        }
//        for(i=0; i < N; i++)
//          Wt[i] /= sum;

//        for(j=0; j < p; j++) {
//          for(k=0; k < p; k++) {
//                        sampcov[j + k*p] = 0;
//                for(i=0; i < N; i++) {
//                sampcov[j + k*p] += Wt[i]*(x[i+N*j]-mu[g+G*j])*(x[i+N*k]-mu[g+G*k]);
//                }
//             }
//           }

//        free(Wt);
//}
void Covariance1(int N, int p, double *x, double *z, double *mu0, double *sampcov) {
  int j,k,i;
        double sum;
        double *Wt = (double*)malloc(sizeof(double)*N);

  sum=0;
        for(i=0; i < N; i++) {
                Wt[i] = z[i];
                sum  += Wt[i];
        }
        for(i=0; i < N; i++)
          Wt[i] /= sum;

        for(j=0; j < p; j++) {
          for(k=0; k < p; k++) {
                        sampcov[j + k*p] = 0;
                for(i=0; i < N; i++) {
                sampcov[j + k*p] += Wt[i]*(x[i+N*j]-mu0[j])*(x[i+N*k]-mu0[k]);
                }
             }
           }

        free(Wt);
}
void Covariance(int N, int p, double *x, double *z, double *mu0, double *sampcov) {
  int j,k,i;
        double sum,  wsum, wsum2;
        double *Wt = (double*)malloc(sizeof(double)*N);

  sum=0;
        for(i=0; i < N; i++) {
                Wt[i]=z[i];
                sum  += Wt[i];
        }
        for(i=0; i < N; i++){
          Wt[i] /= sum;
        }

        for(j=0; j < p; j++) {
          for(k=0; k < p; k++) {
                        wsum = 0.0;
                        wsum2 = 0.0;
                        sampcov[j + k*p] = 0;
                for(i=0; i < N; i++) {
                sampcov[j + k*p] += Wt[i]*(x[j+i*p]-mu0[j])*(x[k+i*p]-mu0[k]);
                    wsum += Wt[i];
                    wsum2 += Wt[i]*Wt[i];
                 }
                   }
                }

        free(Wt);
}
void svd(int M, int N, double *A, double *s, double *u, double *vtt) {
#define LDA M               //M is number of rows of A
#define LDU M
#define LDVT N              //N is number of columns of A
        int i, j;
        int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info, lwork;
        double wkopt;
        double* work;
        double vt[LDVT*N];
        lwork = -1;
        dgesvd_( "All", "All", &m, &n, A, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork,
         &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        dgesvd_( "All", "All", &m, &n, A, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
         &info );
//      for(i=0;i<N;i++)
//          for(j=0;j<N;j++)
//              vtt[i*N+j] = vt[i+j*N];
      for(i=0;i<N;i++)
          for(j=0;j<N;j++)
              vtt[i*N+j] = vt[i*N+j];

        free( (void*)work );
}

void svd1(int M, int N, double *A, double *s, double *u, double *vtt) {
#define LDA M               //M is number of rows of A
#define LDU M
#define LDVT N    

        int i, j;
        int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info, lwork;
        double wkopt;
        double* work;
        int iwork[8*N];
        double vt[LDVT*N];
        lwork = -1;
        dgesdd_( "Singular vectors", &m, &n, A, &lda, s, u, &ldu, vt, &ldvt, &wkopt,
        &lwork, iwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        dgesdd_( "Singular vectors", &m, &n, A, &lda, s, u, &ldu, vt, &ldvt, work,
        &lwork, iwork, &info );
      for(i=0;i<N;i++)
          for(j=0;j<N;j++)
              vtt[i+N*j] = vt[i*N+j];

       free( (void*)work );
}


void eigen(int N, double *A, double *wr, double *vr) {
    int lda = N, ldvl = N, ldvr = N, info, lwork;
    double wkopt;
    double* work;
    double  *wi = malloc(sizeof(double)*N);
     double vl[ldvl*N];

    /* Query and allocate the optimal workspace */
     lwork = -1;
     dgeev_( "Vectors", "Vectors", &N, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     &wkopt, &lwork, &info );
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     /* Solve eigenproblem */
     dgeev_( "Vectors", "Vectors", &N, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     work, &lwork, &info );
      /* Check for convergence */
        if( info > 0 ) {
                Rprintf( "The algorithm failed to compute eigenvalues.\n" );
//               exit( 1 );
                       }
                       //        print_eigenvalues( "Eigenvalues", N, wr, wi );
                      //        print_eigenvectors( "Left eigenvectors", N, wi, vl, ldvl );
                       //        print_eigenvectors( "Right eigenvectors", N, wi, vr, ldvr );
                               dgeev_sort(wr, wi, vr, N);
                                       free(wi);
                                               free( (void*)work );
                                               //        exit( 0 );
                                               }
//

void dgeev_sort(double *Er, double *Ei, double *vr, int N)
{
        double temp;
        double *E2 = malloc(sizeof(double)*N);
        int i, j, k;

        for (i=0; i<N; i++)
                E2[i] = Er[i]*Er[i]+Ei[i]*Ei[i];

        for (j=0; j<N; j++)
                for (i=0; i<N-1; i++)
                        if (fabs(E2[i])<fabs(E2[i+1]))
                        {
                                temp = E2[i]; E2[i] = E2[i+1]; E2[i+1] = temp;
                                temp = Er[i]; Er[i] = Er[i+1]; Er[i+1] = temp;
                                temp = Ei[i]; Ei[i] = Ei[i+1]; Ei[i+1] = temp;

                                for (k=0; k<N; k++)
                                {
                                        temp = vr[k + i*N];
                                        vr[k + i*N] = vr[k + (i+1)*N];
                                        vr[k + (i+1)*N] = temp;
                                }
                        }

       free(E2);
}


// compute the pseudo-inverse of the matrix A and store it in B
int ginv(__CLPK_integer n, __CLPK_integer lda, double *A, double *B) {
        int i;
        __CLPK_integer info;
        __CLPK_integer NRHS = n;
        __CLPK_integer LDWORK = -1;
        __CLPK_integer LDB= n;
        __CLPK_integer IRANK = -1;
        double RCOND = -1.0f;
        double *WORK = (double *)malloc(sizeof(double));
        double *sing = (double *)malloc(sizeof(double)*n);
        double *C = (double *) malloc(sizeof(double)*n*n);

        for(i=0; i < n*n; i++) {
                if(i/n == i%n )
                        B[i] = 1.0;
                else
                        B[i] = 0.0;
        }

        copymx(A, n, n, lda, C);

        dgelss_(&n, &n, &NRHS, C, &n, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &info);
        if(COMMENTS && info != 0) {
                Rprintf("Failed in computing pseudo inverse.\n");
        } else {
                LDWORK = WORK[0];
          free(WORK);
          WORK = (double *)malloc(sizeof(double)*LDWORK);
          dgelss_(&n, &n, &NRHS, C, &n, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &info);

          if(COMMENTS && info != 0) {
                  Rprintf("Failed in computing pseudo inverse.\n");
          }
        }
        free(WORK);
        free(sing);
        free(C);

        return info;
}

void copymx(double *A, int r, int c, int lda, double *C) {
  int i,j;
  for(i=0; i < r; i++) {
        for(j=0; j < c; j++)
                C[i + j*r] = A[i + j*lda];
  }
}


void variance(double *x, int N, int p, double *var){
     int i,j;
     double *mean= malloc(sizeof(double)*p);
     for(j=0; j<p; i++){
        mean[j] =0.0;
        for(i=0; i<N; j++){
           mean[j] += x[j+i*p];
        }
        mean[j] /=N;
     }
     for(j=0; j<p; i++){
        for(i=0; i<N; j++){
           var[j] +=pow((x[j+i*p]-mean[j]),2)/((double) N - (double)1.0);
        }
        var[j]=sqrt(var[j]);
     }
     free(mean);
}
/*void mahalanobis( int N, int p, double *x, double *mu0, double *cov, double *delta0)
{
    int i, j, n;
    double sum, insum;
    double *inv = (double *) malloc(sizeof(double)*p*p);
    ginv(p, p, cov, inv);

    for (n = 0; n < N; n++){
              sum = 0.0;
         for(j = 0; j < p; j++){
            insum = 0;
               for(i = 0; i < p; i++)
                  insum += (x[n+ N*i] - mu0[i])*cov[i+j*p];
//                  insum += (x[n+ N*i] - mu0[i])*inv[i+j*p];
                   sum += insum*(x[n+ N*j] - mu0[j]);
            }
       delta0[n] = sum;
       }
       free(inv);
}*/

     void mahalanobis(int N, int p, double *x, double *mu0, double *cov, double *delta){

      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double *xmu    = malloc(sizeof(double)*N*p);
      double *xmucov    = malloc(sizeof(double)*N*p);
      double *xmucovx    = malloc(sizeof(double)*N*p);
      for (i=0;i<N;i++)
          for(j=0;j<p;j++)
//              xmu[i*p+j]=x[i+j*N]-mu0[j];
              xmu[i*p+j]=x[i*p+j]-mu0[j];

    double *inv = (double *) malloc(sizeof(double)*p*p);
    ginv(p, p, cov, inv);
//    dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,inv,&p,xmu,&p,&beta,xmucov,&p);
    dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,cov,&p,xmu,&p,&beta,xmucov,&p);
      for (i=0;i<N;i++)
          for(j=0;j<p;j++)
             xmucovx[i*p+j]=xmucov[i*p+j]*xmu[i*p+j];
      for (i=0;i<N;i++){
          delta[i]=0.0;
          for(j=0;j<p;j++){
             delta[i] +=xmucovx[i*p+j];
          }
      }
      free(xmu); free(xmucov); free(xmucovx);
}



int determinant(double *A, __CLPK_integer k, __CLPK_integer lda, double *res) {
        int i;
        __CLPK_integer *ipiv;
        double *C;
        __CLPK_integer info=0;

        C = (double*)malloc(sizeof(double)*k*k);
        copymx(A, k, k, lda, C);

        ipiv = (__CLPK_integer *)malloc(sizeof(__CLPK_integer)*k);

        dgetrf_(&k, &k, C, &k, ipiv, &info);

        if(COMMENTS && info != 0)
                Rprintf("Failed in computing matrix determinant.\n");

        *res=1.0;
        for(i=0; i < k*k; i++) {
          if(i%k == i/k)
                  (*res) *= C[i];
        }

        if(*res < 0)     //computing absolute value of the determinant
          *res = -(*res);

        free(ipiv);
        free(C);

        return info;
}


           
/*void iparM(double *x, int N, int p, double *wt, int q, int G,double *mu0, double *alpha, double *Lambda, double *sigma, double *err, double *cpl){
      int i,g,j;
//      char trans = 'T';
      char notrans = 'N';
      char trans = 'T';
      double gamma = 1.0f;
      double beta = 0.0f;
      double *vr= malloc(sizeof(double)*p*p);
      double *wr= malloc(sizeof(double)*p);
      double *var= malloc(sizeof(double)*p);
      double *wmean = malloc(sizeof(double)*p);
      double *dummy    = malloc(sizeof(double)*p*q);
      double *Lambdat    = malloc(sizeof(double)*q*p);
      double **dia    = malloc(sizeof(double*)*G);
      double **sampcov    = malloc(sizeof(double*)*G);
      for(g=0; g<G; g++){
         dia[g]      = malloc(sizeof(double)*q*q);
         sampcov[g]      = malloc(sizeof(double)*p*p);
      }
      if(!wt){
        for(i=0; i<N; i++)
           wt[i] = (double)1.0;
      }    
      weightedmean(x,N,p,wt,wmean);
         for(j=0; j<p; j++)
         mu0[j] = randn(wmean[j],sqrt((double)1.0/(double)N));
      for(g=0; g<G; g++)
         for(j=0; j<p; j++)
         alpha[g*p+j] = randn(0.0,sqrt((double)1.0/(double)N));
         Covariance1(N, p, x, wt, mu0, sampcov[g]);
         eigen(p,sampcov[g],wr,vr);
         for(i=0; i<p*q; i++){
            Lambda[i]=vr[i];
         }
      
         for(i=0; i<q*q; i++)
            dia[g][i]=0.0;
         Covariance1(N, p, x, wt, mu0, sampcov[g]);
         eigen(p,sampcov[g],wr,vr);
         for(i=0; i<q; i++){
            dia[g][i*q+i]=wr[i];
         }
     
      if(q==1){
         dgemm_(&notrans,&trans,&p,&p,&q,&gamma,Lambda,&p,Lambda,&q,&beta,sigma,&p);
        }
      else{
          dgemm_(&notrans,&notrans,&q,&p,&q,&gamma,dia[g],&q,Lambda,&q,&beta,dummy,&q);
        for(i=0;i<q;i++)
           for(j=0;j<p;j++)
              Lambdat[i*p+j]=Lambda[i+j*q];
       dgemm_(&notrans,&notrans,&p,&p,&q,&gamma,Lambdat,&p,dummy,&q,&beta,sigma,&p);
}
     
         for(i=0; i<p*p; i++)
            err[i]= 0.0;
          
         Covariance1(N, p, x, wt, mu0, sampcov[g]);
         for(i=0; i<p; i++){
            err[i*p+i]= sampcov[g][i*p+i]-sigma[i*p+i];
         }
      
          
          
         eigen(p,sigma,wr,vr);
         for(i=0; i<p; i++){
            if(wr[i]<=0.0){
              for(j=0; j<p*p; j++){
              sigma[j] =0.0;
             }
             break;
           }
        }
              variance(x,N,p,var);
              for(j=0; j<p; j++)
                 sigma[j*p+j]=var[j];
             
          
       
//      for(i=0; i<2; i++)
         for(g=0; g<G; g++){
            cpl[0+g*2]= 1.0;
            cpl[1+g*2]= -0.5;
         }
}*/

/*      void rmpar (int p, double *l, int G){
      int i,j,g;
      double cpl0[2];
      double wg1;
      double wg[2]; 
      double *vr= malloc(sizeof(double)*p*p);
      double *cpl= malloc(sizeof(double)*p*2);
      double *wr= malloc(sizeof(double)*p);
      double *mu    = malloc(sizeof(double)*G*p);
      double *phi    = malloc(sizeof(double)*G*p);
      double *alpha    = malloc(sizeof(double)*G*p);
//      double *cpl = malloc(sizeof(double)*p*2);
      double **gam    = malloc(sizeof(double*)*G);
      double **sigma    = malloc(sizeof(double*)*G);
      double **matrix    = malloc(sizeof(double*)*G);
      double **cov1    = malloc(sizeof(double*)*G);
      for(g=0; g<G; g++){
         gam[g]      = malloc(sizeof(double)*p*p);
         cov1[g]      = malloc(sizeof(double)*p*p);
         sigma[g]      = malloc(sizeof(double)*p*p);
         matrix[g]      = malloc(sizeof(double)*(p+1)*p);
      }
      for(g=0; g<G; g++)
         for(j=0; j<p; j++)
            mu[g*p+j]=l[g*p+j];

      for(g=0; g<G; g++)
         for(j=0; j<p; j++)
            phi[g*p+j]=1.0;

      for(g=0; g<G; g++)
         for(j=0; j<p; j++)
         alpha[g*p+j] = randn(0.0,.01);

      for(i=0; i<p; i++){
          cpl[0+i*2]=1.0;
          cpl[1+i*2]=-0.5;
      }

      for(g=0; g<G; g++)
         for(j=0; i<p+1; j++)
             for(j=0; j<p; j++)
                matrix[g][i*p+j]=randn(0,1);

      for(g=0; g<G; g++){
      cov(matrix[g],cov1[g],p+1,p);
         eigen(p,cov1[g],wr,vr);
         for(i=0; i<p*p; i++){
            gam[g][i]=vr[i];
         }
       }

      for(g=0; g<G; g++)
         for(i=0; i<p*p; i++)
             sigma[g][i]=0.0;

      for(g=0; g<G; g++)
         for(i=0; i<p; i++)
             sigma[g][i*p+i]=1.0;
 
      cpl0[0] = 1.0;
      cpl0[1] = -0.5;
      time_t t; 
      srand((unsigned) time(&t));
      wg1 = (double) rand () / RAND_MAX;
      wg[0] = wg1;
      wg[1] = 1.0-wg1;
//      for(g=0; g<G; g++){
//         free(gam[g]);
//         free(cov1[g]);
//         free(sigma[g]);
//         free(matrix[g]);
//      }
//      free(gam); free(cov1); free(sigma); free(matrix);
//      free(vr); free(wr); free(mu); free(phi); free(alpha);
//      free(cpl); 
}*/

      void weights(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double *wg, double *phi, int N, int G, int p,int v, double *pi, double *w){  
      int g,i,j;
      double cpl00[2];
      double wg0[2];
      int n=2;
      double *zlog= malloc(sizeof(double)*N*G);
      double *wt= malloc(sizeof(double)*G*N);
      double *mu0 = malloc(sizeof(double)*p);
      double *alpha0 = malloc(sizeof(double)*p);
      double *phi0 = malloc(sizeof(double)*p);
      double *zlog0 = malloc(sizeof(double)*N);
      double sum =0.0; 
      for(g=0; g<G; g++){
          for(j=0; j<p; j++){
              mu0[j]=mu[g*p+j];
              alpha0[j] = alpha[g*p+j];
              phi0[j] = phi[g*p+j];
          }
          for(i=0;i<n;i++){
             cpl00[i]= cpl0[g*n+i];
             wg0[i] = wg[g*n+i];
           }
         ddmsghyp(x, mu0, alpha0,phi0, cpl00, cpl[g], gam[g], wg0, zlog0, N, G, p);
         for(j=0;j<N;j++){
             zlog[g+j*G] = zlog0[j];
         }
     }
     if(G>1){
        for(i=0;i<N;i++){
            sum=0.0;
            for(g=0;g<G;g++){
               wt[i+g*N]=pow((zlog[i*G+g]*pi[g]),v);
               sum +=wt[i+g*N];
             }
          for(g=0;g<G;g++)
             wt[i+g*N] /=sum;
        if(sum==0){
          for(g=0;g<G;g++)
             wt[i+g*N]=(double)1.0/G;
        }
     }
     for(i=0;i<N;i++)
        for(g=0;g<G;g++)
           w[i*G+g]=wt[i+g*N];
    }
     else{
         for(i=0;i<N;i++)
            for(g=0;g<G;g++)
               w[i*G+g]=(double)1.0;
          }
       

       

     free(zlog); free(wt); free(mu0); free(alpha0); free(zlog0); free(phi0);
}      

      void ddmsghyp(double *x, double *mu0, double *alpha0, double *phi0, double *cpl0, double *cpl, double *gam, double *wg, double *zlog0, int N, int G, int p){

     int i;
     double *val1= malloc(sizeof(double)*N);
     double *val2= malloc(sizeof(double)*N);
         ddghyp(x, mu0, alpha0, cpl0, gam, val1, N, p,phi0);
         dmsghyp(x, mu0, alpha0, phi0, cpl, gam, val2,N,p);
     for(i=0;i<N;i++){
        zlog0[i]=wg[0]*exp(val1[i])+wg[1]*exp(val2[i]);
     }
     free(val1); free(val2);
}
      void ddghyp(double *y, double *mu0, double *alpha0, double *cpl0, double *gam,  double *val1, int N, int p, double *phi0){
     int i,j;
     int m=3;
     double omega, lambda;
     double pa;
     char notrans = 'N';
     double gamma = 1.0f;
     double beta = 0.0f;
     double sum=0;
     double dummy2=0.0;
     double lv[5];
     double *x    = malloc(sizeof(double)*N*p);
     double *xmu  = malloc(sizeof(double)*N*p);
     double *dummy3  = malloc(sizeof(double)*N*p);
     double *dummy1  = malloc(sizeof(double)*p);
     double *dummy4  = malloc(sizeof(double)*N);
     double *lvx  = malloc(sizeof(double)*N*m);
     double *delta0= malloc(sizeof(double)*N);
     double *xmuphi= malloc(sizeof(double)*N);
     double *mx= malloc(sizeof(double)*N);
     double *kx= malloc(sizeof(double)*N);
     double *invS = malloc(sizeof(double)*p*p);
     dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,gam,&p,y,&p,&beta,x,&p);
     omega = cpl0[0];
     lambda= cpl0[1];
//     ginv(p,p,sigma,invS);
//     vec_mx_mult(p,p,alpha0,invS,dummy1);
//      for(i=0;i<p;i++)
//      dummy2  +=dummy1[i]*alpha0[i];
     dummy2 =0.0;
     for(i=0;i<p;i++)
        dummy2 +=pow(alpha0[i],2)/phi0[i];
      pa = omega + dummy2;
      for(i=0; i<N; i++)
         for(j=0; j<p; j++)
            xmu[i*p+j] = x[i*p+j]-mu0[j];
      for(i=0; i<N; i++){
         xmuphi[i]=0.0;
         for(j=0; j<p; j++){
            xmuphi[i] += pow(xmu[i*p+j],2)*1.0/phi0[j];
         }
      }
//      mahalanobis(N,p,x,mu0,invS,delta0);
      for(i=0; i<N; i++){
          mx[i]= omega+xmuphi[i];
          kx[i]=sqrt(mx[i]*pa);
      }

//      dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,invS,&p,xmu,&p,&beta,dummy3,&p); 
//      mx_vec_mult(N,p,alpha0,dummy3,dummy4); 
      for(i=0; i<N; i++){
         dummy4[i]=0.0;
         for(j=0; j<p; j++){
             dummy4[i] +=xmu[i*p+j]*alpha0[j]/phi0[j];
         }
      }
      double nu=lambda -(double)p/(double)2.0;
      int exponscaled = 1;
      int nz = 1;
      int nSeq =1;
      double *r0= malloc(sizeof(double)*nz);
      double *r= malloc(sizeof(double)*N);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
        BesselK(kx[i],nz,nu,exponscaled,nSeq,r0,r1);
        r[i]=r0[0];
      }
      for(i=0; i<N; i++){
         lvx[0+i*m]=(lambda-((double)p/(double)2.0))*log(kx[i]);     
         lvx[1+i*m]= log(r[i])-kx[i];    
         lvx[2+i*m]= dummy4[i];    
      }
     nz =1;
     nu=lambda;
     double *r2= malloc(sizeof(double)*nz);
     BesselK(omega,nz,nu,exponscaled,nSeq,r2,r1);
     double sumphi=0.0;
     for(j=0; j<p; j++)
        sumphi +=log(phi0[j]);
//     determinant(sigma, p, p, &det);
     lv[0]=-((double)1.0/(double)2.0)*sumphi-((double)p/(double)2.0)*(log(2.0)+log(M_PI));

     lv[1]=omega-log(r2[0]);
     lv[2]=-(lambda/(double)2.0)*log(1);
     lv[3]=lambda*log(omega)*0;
     lv[4]=((double)p/(double)2-lambda)*log(pa); 
     for(i=0;i<5;i++)
        sum = sum+lv[i];
     for(i=0;i<N;i++){
        val1[i]=0;
        for(j=0;j<m;j++){
           val1[i] +=lvx[i*m+j];
        }
     val1[i]=val1[i]+sum;
     }
     
     free(x); free(xmu); free(dummy3); free(dummy1); free(dummy4);
     free(lvx); free(delta0); free(mx); free(kx); free(invS); 
     free(r0); free(r); free(r1); free(r2); free(xmuphi);
}
           
      void dmsghyp(double *y, double *mu0, double *alpha0, double *phi0, double *cpl, double *gam, double *val2, int N, int p){
     int i,j;
     int nz;     
     char notrans = 'N';
     double gamma = 1.0f;
     double beta = 0.0f;
     double *x    = malloc(sizeof(double)*N*p);
     double *m0    = malloc(sizeof(double)*N*p);
     double *mx    = malloc(sizeof(double)*N*p);
     double *kx    = malloc(sizeof(double)*N*p);
     double *lx1    = malloc(sizeof(double)*N*p);
     double *lx2    = malloc(sizeof(double)*N*p);
     double *matrix    = malloc(sizeof(double)*N*p);
     double *sum    = malloc(sizeof(double)*N);
     double *lvsum    = malloc(sizeof(double)*p);
     double *lx3    = malloc(sizeof(double)*N*p);
     double *lv1    = malloc(sizeof(double)*p);
     double *lv2   = malloc(sizeof(double)*p);
     double *lv3    = malloc(sizeof(double)*p);
     double *pa    = malloc(sizeof(double)*p);
     double *row    = malloc(sizeof(double)*p);
     double *chi    = malloc(sizeof(double)*p);
     double *psi    = malloc(sizeof(double)*p);
     double *lambda    = malloc(sizeof(double)*p);
     double *dummy  = malloc(sizeof(double)*N*p);
     double *xmu  = malloc(sizeof(double)*N*p);
      double *nu1= malloc(sizeof(double)*p);
     dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,gam,&p,y,&p,&beta,x,&p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            xmu[i*p+j] = x[i*p+j]-mu0[j];
     for(i=0;i<p;i++){
        chi[i]=cpl[0+i*2];
        psi[i]=cpl[0+i*2];
        lambda[i]=cpl[1+i*2];
     }
     for(i=0;i<p;i++)
        pa[i]=psi[i]+pow(alpha0[i],2.0)/phi0[i];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            m0[i*p+j]=pow(xmu[i*p+j],2.0)*(double)1.0/phi0[j];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           mx[i*p+j]=m0[i*p+j]+chi[j];   
     for(i=0;i<N;i++){
        for(j=0;j<p;j++){
           kx[i*p+j]=sqrt(mx[i*p+j]*pa[j]);
         }
     }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           dummy[i*p+j]=log(mx[i*p+j])-log(pa[j]);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           lx1[i*p+j]=dummy[i*p+j]*(lambda[j]-(double)1.0/(double)2.0)/(double)2.0;
     int exponscaled=1;
     nz = 1;
     int nSeq = 1;
     for(j=0;j<p;j++)
         nu1[j] = lambda[j] -(double)1/(double)2;
      double *r0= malloc(sizeof(double)*nz);
      double *r= malloc(sizeof(double)*p);
      double *r1= malloc(sizeof(double)*nz*nSeq);

     for(i=0;i<N;i++){
        for(j=0;j<p;j++){
           row[j]=kx[i*p+j];
           BesselK(row[j],nz,nu1[j],exponscaled,nSeq,r0,r1);
           r[j]=r0[0];
         }
        for(j=0;j<p;j++){
           lx2[i*p+j]=log(r[j])-kx[i*p+j];
        }
     }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            lx3[i*p+j]=xmu[i*p+j]*alpha0[j]/phi0[j];
     for(i=0;i<p;i++)
        lv1[i]=-(double).5*log(phi0[i])-(double).5*(log((double)2)+log(M_PI));
        for(j=0;j<p;j++){
           BesselK(chi[j],nz,lambda[j],exponscaled,nSeq,r0,r1);
           r[j]=r0[0];
        }
     for(i=0;i<p;i++){
        lv2[i]=-(log(r[i])-chi[i]);
         lv3[i]=0;
     }
     for(i=0;i<N;i++) 
        for(j=0;j<p;j++)
           matrix[i*p+j]=lx1[i*p+j]+lx2[i*p+j]+lx3[i*p+j];
     for(i=0;i<N;i++){
         sum[i]=0.0;
         for(j=0;j<p;j++){
            sum[i] +=matrix[i*p+j];
         }
      }
      double lvtot=0.0;
      for(j=0;j<p;j++){
         lvsum[j]=lv1[j]+lv2[j]+lv3[j];
         lvtot +=lvsum[j];
      }
     for(i=0;i<N;i++){
        val2[i]= sum[i]+lvtot;
     }
      
     free(x); free(m0); free(mx); free(kx); free(lx1); free(lx2);
     free(matrix); free(sum); free(lvsum); free(lx3); free(lv1);
     free(lv2); free(lv3); free(pa); free(row); free(chi); free(psi);
     free(lambda); free(dummy); free(xmu); free(nu1); free(r0);
     free(r); free(r1);

}

      void iweights(double *x, double *mu0, double *alpha0, double *phi0, double *cpl, double *cpl0, double *gam, double *wg, int N, int p,int v, double *vg){
      int i,g;
      double sum;
      int G=2;
      double *zlog= malloc(sizeof(double)*N*G);
      double *wt= malloc(sizeof(double)*G*N);
      double *zlog1= malloc(sizeof(double)*N);
      double *zlog2= malloc(sizeof(double)*N);
      ddghyp(x, mu0, alpha0, cpl0, gam, zlog1,N,p,phi0);
      dmsghyp(x, mu0, alpha0, phi0, cpl, gam, zlog2,N,p);

      for(i=0;i<N;i++){
         zlog1[i]=exp(zlog1[i]);
         zlog2[i]=exp(zlog2[i]);
      }
      for(i=0;i<N;i++){
          zlog[0+i*G]=zlog1[i];
          zlog[1+i*G]=zlog2[i];
      }
     if(wg[1] != 0){
     for(i=0;i<N;i++){
        sum=0.0;
        for(g=0;g<G;g++){
            vg[i*G+g]=pow(zlog[i*G+g]*wg[g],v);
           sum +=vg[i*G+g];
        }
        for(g=0;g<G;g++){
           vg[i*G+g] = vg[i*G+g]/sum;
        }
        if(sum==0){
          for(g=0;g<G;g++)
             vg[i*G+g]=(double)1.0/G;
        }
     }
     }else{
     for(i=0;i<N;i++){
        vg[0+i*2]=1.0;
        vg[1+i*2]=0.0;
    }
   }
         
          
     free(zlog); free(wt); free(zlog1); free(zlog2); 
}


void combinewk (double *z, int N, int G, int *labels)
{
     int i, g, n;
     int sum =0;
     for(i = 0; i < N; i++)
     {  if(labels[i] ==0)
           {}
        else
        sum = sum +1;
     }
     for(n=0; n<sum; n++)
         for(g=0; g<G; g++){
             z[n + N*g] = 0.0;
         }
       for (i = 0; i<N; i++){
           if(labels[i] ==0){}
           else
               z[i + N*(labels[i]-1)] =  (double)1.00;
           }
       
}


      void EMgrstep(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double *wg, double *phi, int N, int G, int p, int v, double *pi, int *label,int m){

      double wg0[2];
      double meanvg[2];
      double cpl00[2];
      double *mu0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *phi0= malloc(sizeof(double)*p);
      double *ng= malloc(sizeof(double)*G);
      double *wt= malloc(sizeof(double)*N);
      double *w= malloc(sizeof(double)*N*G);
      double *vg= malloc(sizeof(double)*N*2);
      int i,j,k,g;
      int n =2;
//      for(g=0;g<G;g++){
         weights(x,mu,alpha,cpl,cpl0,gam,wg,phi,N,G,p,v,pi,w); 
      if(label)
        combinewk1(w,N,G,label);
      for(k=0;k<G;k++){
         ng[k]=0.0;
         for(i=0;i<N;i++){
            ng[k] +=w[k+i*G];
         }
      }

      for(g=0;g<G;g++){
         for(j=0;j<p;j++){
            mu0[j]=mu[g*p+j];
            alpha0[j]=alpha[g*p+j];
            phi0[j]=phi[g*p+j];
         }
         for(i=0;i<n;i++){
            cpl00[i]= cpl0[g*n+i];
            wg0[i]= wg[g*n+i];
        }
         iweights(x,mu0,alpha0,phi0,cpl[g],cpl00,gam[g],wg0,N,p,v,vg);
      for(i=0;i<N;i++){
         wt[i]=w[g+i*G];}
         
      updatemaScplp(x,mu0,gam[g],phi0,alpha0,N,p,wt,cpl[g],cpl00,vg,m);
      weightedmean(vg,N,2,wt,meanvg);
      if(meanvg[1] < 0.01){
         meanvg[0]=.99;
         meanvg[1]=.01;
        }
         for(i=0;i<n;i++){
            wg[g*n+i] = meanvg[i];
         }
      for(j=0;j<p;j++){
         mu[g*p+j]= mu0[j];
         alpha[g*p+j] = alpha0[j];
         phi[g*p+j]= phi0[j];
      }
      for(i=0;i<n;i++){
         cpl0[g*n+i]=cpl00[i];
      } 
      }
      for(g=0;g<G;g++){
         pi[g]=0.0;
         for(i=0;i<N;i++){
            pi[g] +=w[g+i*G];
         }
         pi[g] /= N;
      }
      free(mu0); free(alpha0); free(phi0); free(ng); free(wt);
      free(w); free(vg);   


}
      void grad2(double *x, int N, double *y, double eps, double d, int r, int v, double *vec){

     int i,j,k;
     int  besselkAsym =0;
     double zerotol=sqrt(2.220446e-16/7e-7);
     double *h= malloc(sizeof(double)*N);
     double *a= malloc(sizeof(double)*r*N);
     for(i=0;i<N;i++){
        h[i]=fabs(d*x[i]);
           if(fabs(x[i]) < zerotol)
             h[i] +=eps;
     }
     int exponscaled=0;
     int nz = 1;
     int nSeq = 1;
     double dummy;
     double res0, res1, res2;
     double *r0= malloc(sizeof(double)*nz);
     double *r1= malloc(sizeof(double)*nz*nSeq);
     for(j=0;j<r;j++){
        for(i=0;i<N;i++){
           dummy=y[i];
           BesselK(dummy,nz,x[i],exponscaled,nSeq,r0,r1);
//           res0= log(r0[0])-y[i];
           res0= log(r0[0]);
           dummy=y[i];
           if(res0 == INFINITY || res0 == -INFINITY)
//             {besselkAsym =1;
//              break;}
           {
           int kmax = 4;
           x[i]=fabs(x[i]);
           besselKnuAsym(dummy,x[i]+h[i],kmax,&res1);
           besselKnuAsym(dummy,x[i]-h[i],kmax,&res2);
           a[j*N+i] = (res1-res2)/(2*h[i]);
           }else{
           dummy=y[i];
           BesselK(dummy,nz,x[i]+h[i],exponscaled,nSeq,r0,r1);
//           res1 = log(r0[0]) -(y[i]);
           res1 = log(r0[0]) ;
           dummy=y[i];
           BesselK(dummy,nz,x[i]-h[i],exponscaled,nSeq,r0,r1);
//           res2 = log(r0[0])-(y[i]);
           res2 = log(r0[0]);
          a[j*N+i] = (res1-res2)/(2*h[i]);
        }
     }      
        for(k=0;k<N;k++){
           h[k] /=v;
        }
//     }
//     for(i=0;i<N;i++){
//        h[i]=fabs(d*x[i]);
//           if(fabs(x[i]) < zerotol)
//             h[i] +=eps;
//     }

//     int kmax=4;
//     if(besselkAsym){
//     for(j=0;j<r;j++){
//        for(i=0;i<N;i++){
//           dummy=y[i];
//           x[i]=fabs(x[i]);
//           besselKnuAsym(dummy,x[i]+h[i],kmax,&res1);
//           besselKnuAsym(dummy,x[i]-h[i],kmax,&res2);
//           a[j*N+i] = (res1-res2)/(2*h[i]);
//        }
//        for(k=0;k<N;k++){
//           h[k] /=v;
//        }
//     }
    }
     for(k=1;k<r;k++){
        double *mat= malloc(sizeof(double)*(r-k)*N);
        for(j=0;j<r-k;j++){
           for(i=0;i<N;i++){
              mat[j*N+i] = (a[(k+j)*N+i]*pow(4,k)-a[(k+j-1)*N+i])/(pow(4,k)-1);
           }
        }
        for(j=0;j<r-k;j++){
           for(i=0;i<N;i++){
              a[(j+k)*N+i]= mat[j*N+i];
              vec[i]= mat[j*N+i];
           }
        }
        free(mat);      
     }
     free(h); free(a); free(r0); free(r1);
}


      void grad(double *x, int N, double *y, double eps, double d, int r, int v, double *vec){

     int i,j,k;
     double zerotol=sqrt(2.220446e-16/7e-7);
     double *h= malloc(sizeof(double)*N);
     double *a= malloc(sizeof(double)*r*N);
     for(i=0;i<N;i++){
        h[i]=fabs(d*x[i]);
           if(fabs(x[i]) < zerotol)
             h[i] +=eps;
        
     }
     int exponscaled=1;
     int nz = 1;
     int nSeq = 1;
     double dummy;
     double res1, res2;
     double *r0= malloc(sizeof(double)*nz);
     double *r1= malloc(sizeof(double)*nz*nSeq);
     for(j=0;j<r;j++){
        for(i=0;i<N;i++){
           dummy=y[i];
           BesselK(dummy,nz,x[i]+h[i],exponscaled,nSeq,r0,r1);
           res1 = log(r0[0]) -log(y[i]);
           dummy=y[i];
           BesselK(dummy,nz,x[i]-h[i],exponscaled,nSeq,r0,r1);
           res2 = log(r0[0])-log(y[i]);
           a[j*N+i] = (res1-res2)/(2*h[i]);
        }
        for(k=0;k<N;k++){
           h[k] /=v;
        }
     }
     for(k=1;k<r;k++){
        double *mat= malloc(sizeof(double)*(r-k)*N);
        for(j=0;j<r-k;j++){
           for(i=0;i<N;i++){
              mat[j*N+i] = (a[(k+j)*N+i]*pow(4,k)-a[(k+j-1)*N+i])/(pow(4,k)-1);
           }
        }
        for(j=0;j<r-k;j++){
           for(i=0;i<N;i++){
              a[(j+k)*N+i]= mat[j*N+i];
              vec[i]= mat[j*N+i];
           }
        }
        free(mat);

     }
     free(h); free(a); free(r0); free(r1);
}

      void besselKnuAsym( double x, double nu, int kmax, double *res){
      double d=0.0;
      double z = x/nu;
//      Rprintf("z is %f %f\n", x,nu);
      double sz=sqrt(1+pow(z,2.0));
      double t = 1.0/sz;
      double eta = sz+log(z/(1+sz));
//      Rprintf(" eta sz z %f %f %f\n", eta, sz, z);
      if(kmax==0) d=0;
        else{
            double t2 = pow(t,2.0);
            double u1t = (t*(3.0-5.0*t2))/24;
            if(kmax==1){
              d=-u1t/nu;
              } else {
                double u2t = t2*(81+t2*(-462 + t2*385))/1152;
                if(kmax==2)
                  d=(-u1t+u2t/nu)/nu;
                  else{ double u3t = t*t2*(30375+t2*(-369603+t2*(765765 - t2*425425)))/ 414720;
                   if(kmax==3)
                     d=(- u1t + (u2t - u3t/nu)/nu)/nu;
                     else{double u4t = t2*t2*(4465125 +t2*(-94121676 +t2*(349922430 +t2*(-446185740 + t2*185910725))))/39813120;
                     if(kmax==4)
                        d=(- u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;
//                        Rprintf("u1t u2t u3t u4t %f %f %f %f\n", u1t, u2t, u3t, u4t);
                    }
                }
             }
         }
         *res=log1p(d)-nu*eta - (log(sz) - log(M_PI/(2*nu)))/2;
//           *res = (1+d)*exp(-nu*eta)*sqrt(M_PI/(2*nu*sz));
//         Rprintf("d nu eta sz MPI res %f %f %f %f %f %f\n", d, nu, eta, sz,M_PI, *res);
}             
       
      

 


      void updatemaScplp(double *y, double *mu0, double *gam, double *phi0, double *alpha0, int N, int p, double *w, double *cpl, double *cpl0, double *vg, int m){
      int i,j,g;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double parol[2];
      double cplnewU[2];
      double ABC[3];
      double *wt1= malloc(sizeof(double)*N);
      double *wt2= malloc(sizeof(double)*N);
      double *A= malloc(sizeof(double)*p);
      double *B= malloc(sizeof(double)*p);
      double *C= malloc(sizeof(double)*p);
      double *S1= malloc(sizeof(double)*N*p);
      double *S2= malloc(sizeof(double)*N*p);
      double *S= malloc(sizeof(double)*N*p);
      double *alphanew= malloc(sizeof(double)*p);
      double *Numal= malloc(sizeof(double)*p);
      double *abc= malloc(sizeof(double)*N*3);
      double *cplnew= malloc(sizeof(double)*p*2);
      double *x= malloc(sizeof(double)*N*p);
      double *xmu= malloc(sizeof(double)*N*p);
      double *munew= malloc(sizeof(double)*p);
      double *sr2= malloc(sizeof(double)*N*p);
      double *Z= malloc(sizeof(double)*N*p);
      double *W= malloc(sizeof(double)*N*p);
      double *invW= malloc(sizeof(double)*N*p);
      double *logW= malloc(sizeof(double)*N*p);
      double *newgam= malloc(sizeof(double)*p*p);
      double *r= malloc(sizeof(double)*N*p);
      double *Num1sig= malloc(sizeof(double)*p*p);
      double *Num2sig= malloc(sizeof(double)*p*p);
      double *phi0new= malloc(sizeof(double)*p);
      double *omega= malloc(sizeof(double)*p);
      double *sS1= malloc(sizeof(double)*p);
      double *sS2= malloc(sizeof(double)*p);
      double *test1= malloc(sizeof(double)*p*2);
      double *test2= malloc(sizeof(double)*p*3);
      double *row1= malloc(sizeof(double)*2);
      double *row2= malloc(sizeof(double)*3);
      double *sumS= malloc(sizeof(double)*p);
      double *moltb= malloc(sizeof(double)*p);
      double *nummu= malloc(sizeof(double)*p);
      double *nummu1= malloc(sizeof(double)*p);
      double *sumalphanew= malloc(sizeof(double)*p);
      double *wvgW= malloc(sizeof(double)*p);
      double *wvgWalphanew= malloc(sizeof(double)*p);
      double *matrix3sum= malloc(sizeof(double)*p);
      double *molt= malloc(sizeof(double)*N*p);
      double *matrix1= malloc(sizeof(double)*N*p);
      double *matrix2= malloc(sizeof(double)*N*p);
      double *matrix3= malloc(sizeof(double)*N*p);
      double *matrix3r= malloc(sizeof(double)*N*p);
      double *matrix2r= malloc(sizeof(double)*N*p);
      double *matrix2rt= malloc(sizeof(double)*p*N);
      double *matrix2rtr= malloc(sizeof(double)*p*p);
      double *alphanewalphanew= malloc(sizeof(double)*p*p);
      double *alphanewalphanew1= malloc(sizeof(double)*p*p);
      double *alphanewalphanew2= malloc(sizeof(double)*p*p);
      double *alphanewalphanew3= malloc(sizeof(double)*p*p);
      double *alphanewalphanew4= malloc(sizeof(double)*p*p);
      double *alphanewalphanew5= malloc(sizeof(double)*p*p);
      double *matrix4= malloc(sizeof(double)*N*p);
      double *matrix5= malloc(sizeof(double)*N*p);
      double *matrix5r= malloc(sizeof(double)*N*p);
      double *matrix5sum= malloc(sizeof(double)*p);
      double *matrix4Wr= malloc(sizeof(double)*N*p);
      double *matrix4Wrt= malloc(sizeof(double)*p*N);
      double *matrix4Wrtr= malloc(sizeof(double)*p*p);
      double *diagp= malloc(sizeof(double)*p*p);
      double *newgamdiagp= malloc(sizeof(double)*p*p);
      double *newgamt= malloc(sizeof(double)*p*p);
      double *meanS1= malloc(sizeof(double)*p);
      double *meanS2= malloc(sizeof(double)*p);
      double *sums= malloc(sizeof(double)*p);
      double *xS= malloc(sizeof(double)*N*p);
      double *xmuS2= malloc(sizeof(double)*N*p);
      double *phinew= malloc(sizeof(double)*p);
//      double *alphanew= malloc(sizeof(double)*p);
      double *xmumean= malloc(sizeof(double)*p);
      double *xmu2= malloc(sizeof(double)*p);
      double *sr2sum= malloc(sizeof(double)*N);
//      double *phinew= malloc(sizeof(double)*p);
      double *S2meanS2= malloc(sizeof(double)*N*p);
      double *xS2meanS2= malloc(sizeof(double)*N*p);
      
      int G=2;
      if(!w)
        for(i=0;i<N;i++)
           w[i]=1;
      if(!vg){
        for(i=0;i<N;i++)
           for(g=0;g<G;g++)
              vg[i*G+g]=(double)0.5;
      }
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam,&p,y,&p,&beta,x,&p);
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
            Z[i*p+j]=pow(x[i*p+j]-mu0[j],2);
            sr2[i*p+j]=Z[i*p+j]*1/phi0[j];
         }
      }
      for(i=0;i<N;i++){
         sr2sum[i]=0.0;
         for(j=0;j<p;j++){
            sr2sum[i] += sr2[i*p+j];
         }
      }
      gig20(sr2sum,N,p,mu0,alpha0,cpl,cpl0,abc,phi0);
      double sumwi[2];  
      for(g=0;g<G;g++){
         sumwi[g]=0.0;
         for(i=0;i<N;i++){
            sumwi[g]+=vg[g+i*G];
         }
         sumwi[g] = sumwi[g]/N;
      }

      double sumw=0.0;
      for(i=0;i<N;i++)
         sumw +=w[i];
      if(sumwi[1] !=0.0){
      gig2p(sr2,N,p,cpl,phi0,alpha0,W,invW,logW);
//      gig2(x,N,p,mu0,alpha0,cpl,cpl0,abc,phi0);
//      time_t t;
//      srand((unsigned) time(&t));
//      double  random = (double) rand () / RAND_MAX;
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             S1[i*p+j]=W[i*p+j]*vg[1+i*2]+abc[0+i*3]*vg[0+i*2];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            S2[i*p+j]=invW[i*p+j]*vg[1+i*2]+abc[1+i*3]*vg[0+i*2];
      }else{
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             S1[j+i*p]=abc[0+i*3];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             S2[j+i*p]=abc[1+i*3];
      }
      weightedmean(S1,N,p,w,meanS1);
      weightedmean(S2,N,p,w,meanS2);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            S[i*p+j] = ((S2[i*p+j]*meanS1[j])-1.0)*w[i];
      for(j=0;j<p;j++){
          sums[j]=0.0;
          for(i=0;i<N;i++){
             sums[j] +=S[j+i*p];
          }
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xS[i*p+j]=x[i*p+j]*S[i*p+j];
      for(j=0;j<p;j++){
          munew[j]=0.0;
          for(i=0;i<N;i++){
             munew[j] +=xS[j+i*p];
          }
          munew[j] = munew[j]/sums[j];
      }
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
            S2meanS2[i*p+j] =-S2[i*p+j]+meanS2[j];
         }
      }

      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xS2meanS2[i*p+j]= x[i*p+j]*S2meanS2[i*p+j];
      for(j=0;j<p;j++){
         alphanew[j]=0.0;
          for(i=0;i<N;i++){
             alphanew[j] +=xS2meanS2[j+i*p]*w[i];
          }
          alphanew[j]= alphanew[j]/sums[j];
     }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             xmu[i*p+j]=x[i*p+j]-munew[j];
      weightedmean(xmu,N,p,w,xmumean);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xmuS2[i*p+j]=pow(xmu[i*p+j],2)*S2[i*p+j];
      weightedmean(xmuS2,N,p,w,xmu2);
     
      for(j=0;j<p;j++){
         phinew[j]=xmu2[j]-2.0*xmumean[j]*alphanew[j]+meanS1[j]*pow(alphanew[j],2);
       }

      for(i=0;i<N;i++)
         wt1[i]=w[i]*vg[1+i*2];
      for(i=0;i<N;i++)
         wt2[i]=w[i]*vg[0+i*2];
      if(sumwi[1] ==0.01){
        weightedmean(W,N,p,wt1,A);
        weightedmean(invW,N,p,wt1,B);
        weightedmean(logW,N,p,wt1,C);

      for(j=0;j<p;j++)
         omega[j]=cpl[0+j*2];
       for(j=0;j<p;j++){
           test1[0+j*2]=omega[j];
           test1[1+j*2]=cpl[1+j*2];
       }
       for(j=0;j<p;j++){
          test2[0+j*3]=A[j];
          test2[1+j*3]=B[j];
          test2[2+j*3]=C[j];
       }
       for(j=0;j<p;j++){
          for(i=0;i<2;i++){
             row1[i]=test1[j*2+i];
          }
          for(i=0;i<3;i++){
             row2[i]=test2[j*3+i];
          }
          updateol(row1,row2,2);
          updateol(row1,row2,2);
          for(i=0;i<2;i++){
             cplnew[j*2+i]=row1[i];
          }
       }
     }else{
       for(j=0;j<p;j++)
          for(i=0;i<2;i++)
             cplnew[j*2+i]=cpl[j*2+i];
     }
     parol[0] = cpl0[0];
     parol[1] = cpl0[1];
      weightedmean(abc,N,3,wt2,ABC);
      updateol(parol,ABC,2);
      updateol(parol,ABC,2);
      cplnewU[0]=parol[0];
      cplnewU[1]=parol[1];
      for(j=0;j<p;j++){
         mu0[j]=munew[j];
         phi0[j]=phinew[j];
         alpha0[j]=alphanew[j];
      }
     if(m%2 ==0){
       updategam1MS(gam, y, phi0, alpha0, mu0, w, S2, newgam, N, p);
      } 
      else{
       updategam2MS(gam, y, phi0, alpha0, mu0, w, S2, newgam, N, p);
      } 
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            gam[i*p+j]=newgam[i*p+j];
      for(i=0;i<p;i++)
          for(j=0;j<2;j++)
             cpl[i*2+j]=cplnew[i*2+j];

       cpl0[0]=cplnewU[0];
       cpl0[1]=cplnewU[1];
       free(A); free(B); free(C); free(S1); free(S2); free(S); free(alphanew);
       free(Numal); free(abc); free(cplnew); free(x); free(xmu); free(munew);
       free(sr2); free(Z); free(W); free(invW); free(logW); free(newgam); 
       free(r); free(Num1sig); free(Num2sig); free(phi0new); free(omega); 
       free(sS1); free(sS2); free(test1); free(test2); free(row1); free(row2);
       free(sumS); free(moltb); free(nummu); free(nummu1); free(sumalphanew);
       free(wvgW); free(wvgWalphanew); free(matrix3sum); free(molt); free(matrix1);
       free(matrix2); free(matrix3); free(matrix3r); free(matrix2r); free(matrix2rt);
       free(matrix2rtr); free(alphanewalphanew); free(alphanewalphanew1); 
       free(alphanewalphanew2); free(alphanewalphanew3); free(alphanewalphanew4);
       free(alphanewalphanew5); free(matrix4); free(matrix5); free(matrix5r);
       free(matrix5sum); free(matrix4Wr); free(matrix4Wrt); free(matrix4Wrtr);
       free(diagp); free(newgamdiagp); free(newgamt);
       free(meanS1); free(meanS2); free(sums); free(xS); free(xmuS2);
       free(phinew); free(xmumean); free(xmu2); free(sr2sum); free(S2meanS2);
       free(xS2meanS2);


}
    


      void updateolU(double *ol, double *ABC, int n){
      double Rp, Rn, f1,f2,bv;
      if(ABC[2]==0){
        ol[1]=0;
      }else{
      double eps = 1e-8;
      double d = 0.0001;
      int r =6;
      int vv=2;
      int N=1;
      double *vec=malloc(sizeof(double)*N);
      double *dummy1=malloc(sizeof(double)*N);
      double *dummy2=malloc(sizeof(double)*N);
      dummy1[0] = ol[1];
      dummy2[0] = ol[0];
        grad(dummy1,N,dummy2,eps,d,r,vv,vec);
      bv = vec[0];
      ol[1]=ABC[2]*ol[1]/bv;
      free(vec); free(dummy1); free(dummy2);
      }
      int nz=1;
      double nu = ol[1];
      int exponscaled = 1;
      int nSeq =1;
      double *r= malloc(sizeof(double)*nz);
      double *r2= malloc(sizeof(double)*nz);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      BesselK(ol[0],nz,nu+1,exponscaled,nSeq,r,r1);
      BesselK(ol[0],nz,nu,exponscaled,nSeq,r2,r1);
      Rp = r[0]/r2[0];
      nu = -ol[1];
       BesselK(ol[0],nz,nu+1,exponscaled,nSeq,r,r1);
       BesselK(ol[0],nz,nu,exponscaled,nSeq,r2,r1);
       Rn=r[0]/r2[0];
       f1 = Rp + Rn - (ABC[0]+ABC[1]);
       f2 = (pow(Rp,2)-(2*ol[1]+1.0)/ol[0]*Rp-1.0)+(pow(Rn,2)-(2.0*(-1*ol[1])+1.0)/ol[0]*Rn-1.0);
       if((ol[0] -f1/f2) > 0)
       ol[0]=ol[0]-f1/f2;
        free(r); free(r2); free(r1);


}


      void gig2p(double *sr2, int N, int p, double *cpl, double *phi0, double *alpha0,double *W, double *invW, double *logW){
      int i,j;
      double *omega= malloc(sizeof(double)*p);
      double *a1= malloc(sizeof(double)*p);
      double *v1= malloc(sizeof(double)*p);
      double *B1= malloc(sizeof(double)*N*p);
      for(i=0;i<p;i++)
         omega[i]=exp(log(cpl[0+i*2]));
      for(i=0;i<p;i++)
         a1[i]=omega[i]+pow(alpha0[i],2)/phi0[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            B1[i*p+j]=sr2[i*p+j]+omega[j];
      for(i=0;i<p;i++)
          v1[i]=cpl[1+i*2]-((double)1.0/(double)2.0);
      gigp(a1,B1,N,p,v1,W,invW,logW);
      free(omega); free(a1); free(v1); free(B1);
}   

      void gigp(double *a, double *B, int N, int p, double *v,double *W, double *invW, double *logW){
      int i,j;
      double *SaB= malloc(sizeof(double)*N*p);
      double *SBa= malloc(sizeof(double)*N*p);
      double *matrix= malloc(sizeof(double)*N*p);
      double *dummy= malloc(sizeof(double)*N*p);
      double *Kv1= malloc(sizeof(double)*N*p);
      double *Kv0= malloc(sizeof(double)*N*p);
      double *Kv12= malloc(sizeof(double)*N*p);
      double *x= malloc(sizeof(double)*N*p);
      double *nu= malloc(sizeof(double)*p);
      double *row= malloc(sizeof(double)*p);
      double *x0= malloc(sizeof(double)*p);
      double *SaB0= malloc(sizeof(double)*p);
      double *vec= malloc(sizeof(double)*p);
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             SaB[j+i*p]=sqrt(B[j+i*p]*a[j]);
      int exponscaled=1;
      int nz =1;
      int nSeq =1;
     for(j=0;j<p;j++)
         nu[j] = v[j]+(double)1 ;
      double *r0= malloc(sizeof(double)*nz);
      double *r= malloc(sizeof(double)*p);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
          for(j=0;j<p;j++){
             row[j]=SaB[i*p+j];
             BesselK(row[j],nz,nu[j],exponscaled,nSeq,r0,r1);
             r[j]=r0[0];
          }
          for(j=0;j<p;j++)
             Kv1[i*p+j]=r[j];
          }
     
     for(j=0;j<p;j++)
         nu[j] = v[j] ;
      for(i=0;i<N;i++){
          for(j=0;j<p;j++){
             row[j]=SaB[i*p+j];
             BesselK(row[j],nz,nu[j],exponscaled,nSeq,r0,r1);
             r[j]=r0[0];
          }
          for(j=0;j<p;j++)
        Kv0[i*p+j]=r[j];
      }
      for(i=0;i<N;i++)
          for(j=0;j<p;j++)
          Kv12[i*p+j]=Kv1[i*p+j]/Kv0[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
     SBa[j+i*p]=sqrt(B[j+i*p]*(double)1/a[j]);
      for(i=0;i<N;i++){
          for(j=0;j<p;j++){
             W[i*p+j]=Kv12[i*p+j]*SBa[i*p+j];
          }
      }
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            dummy[j+i*p]=((double)1.0/B[j+i*p])*v[j];
      for(i=0;i<N;i++)
          for(j=0;j<p;j++)
             invW[i*p+j]=(Kv12[i*p+j]/SBa[i*p+j])-(double)2.0*dummy[i*p+j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            x[i*p+j]=v[j];
      double eps =1e-8;
      double d = 0.0001;
      int rr = 6;
      int vv=2;
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
            x0[j]=x[i*p+j];
            SaB0[j]=SaB[i*p+j];
          }
           grad2(x0,p,SaB0,eps,d,rr,vv,vec);
           for(j=0;j<p;j++){
               matrix[i*p+j]=vec[j];
           }
       }
      for(i=0;i<N;i++)
          for(j=0;j<p;j++)
             logW[i*p+j]=log(SBa[i*p+j])+matrix[i*p+j];
     free(SaB); free(SBa); free(dummy); free(Kv1); free(Kv0);
      free(Kv12); free(x); free(nu); free(row); free(x0); free(SaB0);
      free(vec); free(r0); free(r); free(r1); free(matrix);
}
      void updategam2(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p, double *newgam){ 
      int i,j;
      double sumw = 0.0; 
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double *x= malloc(sizeof(double)*N*p);
      double *wvg1= malloc(sizeof(double)*N);
      double *sign= malloc(sizeof(double)*p);
      double *wvg2= malloc(sizeof(double)*N);
      double *Bs= malloc(sizeof(double)*N*p);
      double *v1= malloc(sizeof(double)*N*p);
      double *v1t= malloc(sizeof(double)*p*N);
      double *v= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *Bs2= malloc(sizeof(double)*N*p);
      double *xBs2= malloc(sizeof(double)*N*p);
      double *xBs2t= malloc(sizeof(double)*p*N);
      double *dummy= malloc(sizeof(double)*N*p);
      double *dummy2= malloc(sizeof(double)*p*p);
      double *dummy3= malloc(sizeof(double)*p*p);
      double *dummy6= malloc(sizeof(double)*p*p);
      double *dummy5= malloc(sizeof(double)*p*p);
      double *s= malloc(sizeof(double)*p);
      double *u= malloc(sizeof(double)*p*p);
      double *ut= malloc(sizeof(double)*p*p);
      double *vtt= malloc(sizeof(double)*p*p);
      double *vttt= malloc(sizeof(double)*p*p);
      double *dummy4= malloc(sizeof(double)*p);
      double *dummyt= malloc(sizeof(double)*p*N);
      double *wty= malloc(sizeof(double)*N*p);
      double *wty2= malloc(sizeof(double)*N*p);
      double *Bsmu= malloc(sizeof(double)*N*p);
      double *e1= malloc(sizeof(double)*N);
      double *e2= malloc(sizeof(double)*N);
      double *F0t= malloc(sizeof(double)*p*p);
      double *F1t= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *A1= malloc(sizeof(double)*p*p);
      double *F2= malloc(sizeof(double)*p*p);
      double *gamt= malloc(sizeof(double)*p*p);
      double *phimt= malloc(sizeof(double)*N*p);
      double *Bs2t= malloc(sizeof(double)*p*N);
      double *Bs2tphimt= malloc(sizeof(double)*p*N);
      double *mumt= malloc(sizeof(double)*N*p);
      double *alphaphi= malloc(sizeof(double)*N*p);
      double *Bs2mualphaphi= malloc(sizeof(double)*p*N);
      double *Bs2mualphaphit= malloc(sizeof(double)*N*p);
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=w[i];
     
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs[i*p+j]=invW[i*p+j]*(double)1.0/phi0[j];
      
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2[j+i*p]=abc[1+i*3];

      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            phimt[i*p+j]=1.0/phi0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2t[j*N+i]=Bs2[j+i*p];
         for(i=0;i<N*p;i++)
            Bs2tphimt[i]=Bs2t[i]*phimt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs2[i*p+j]=Bs2tphimt[i+j*N];
      for(i=0;i<N;i++){
         wvg1[i]=w[i]*vg[0+i*2];
         wvg2[i]=w[i]*vg[1+i*2];
      }
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
            wty[i*p+j]=wvg2[i]*y[i*p+j];
            wty2[i*p+j]=wvg1[i]*y[i*p+j];
         }
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            dummy[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            dummyt[j*N+i]=dummy[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,dummyt,&N,&beta,dummy2,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs2[i*p+j]=x[i*p+j]*Bs2[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            xBs2t[j*N+i]=xBs2[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty2,&p,xBs2t,&N,&beta,dummy3,&p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F0t[i*p+j]=dummy2[i*p+j]+dummy3[i*p+j];
      for(i=0;i<N;i++){
         e1[i]=0.0;
         e2[i]=0.0;
         for(j=0;j<p;j++){
             e1[i] +=pow(y[i*p+j],2)*wvg1[i];
             e2[i] +=pow(y[i*p+j],2)*wvg2[i];
         }
      }
      for(j=0;j<p;j++){
         dummy4[j]=0.0;
         for(i=0;i<N;i++){
            dummy4[j] +=Bs[j+i*p]*e2[i];
         }
      }
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            gamt[i*p+j]=gam0[i+j*p];
      for(i=0;i<p;i++){
         for(j=0;j<p;j++){
            dummy5[i*p+j]=0;
            if(i==j)
              dummy5[i*p+i]=dummy4[i];
         }
     }
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gamt,&p,dummy5,&p,&beta,dummy6,&p);
      double sum=0.0;
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             sum +=Bs2[j+i*p]*e1[i];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1t[i*p+j]=dummy6[i*p+j]+sum*gamt[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            v[j+i*p]=Bs[j+i*p]*mu0[j]+alpha0[j]/phi0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            mumt[i*p+j]=mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             alphaphi[i*p+j]=alpha0[j]/phi0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2t[j*N+i]=Bs2[j+i*p];
      for(i=0;i<N*p;i++)
         Bs2mualphaphi[i]=Bs2t[i]*mumt[i]+alphaphi[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs2mualphaphit[i*p+j]=Bs2mualphaphi[i+j*N];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v1[i*p+j]=Bs2mualphaphit[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            vt[j*N+i]=v[j+i*p];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            v1t[j*N+i]=v1[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,vt,&N,&beta,A0,&p);
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty2,&p,v1t,&N,&beta,A1,&p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F2[i*p+j]=(F0t[i*p+j]-F1t[i*p+j]-A0[i*p+j]-A1[i*p+j])/sumw;
      svd(p,p,F2,s,vtt,u);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            ut[i*p+j]=u[i+j*p];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            vttt[i*p+j]=vtt[i+j*p];
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,ut,&p,vttt,&p,&beta,newgam,&p);
      for(i=0;i<p;i++){
          if(newgam[i*p+i]<0.0)
            sign[i]=-1.0;
          else
             sign[i]=1.0;
       }
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
             newgam[i*p+j]=newgam[i*p+j]*sign[j];
      if(objgam(newgam,y,phi0,alpha0,mu0,w,vg,invW,abc,N,p) <objgam(gam0,y,phi0,alpha0,mu0,w,vg,invW,abc,N,p)){
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
             newgam[i*p+j]=gam0[i*p+j];
       }
      free(x); free(wvg1); free(sign); free(wvg2); free(Bs); free(v);
      free(vt); free(Bs2); free(xBs2); free(xBs2t); free(dummy); free(dummy2);
      free(dummy3); free(dummy6); free(dummy5); free(s); free(u); free(ut);
      free(vtt); free(dummy4); free(dummyt); free(wty); free(wty2); 
      free(Bsmu); free(e1); free(e2); free(F0t); free(F1t); free(A0);
      free(F2); free(gamt); free(phimt); free(Bs2t); free(Bs2tphimt);
      free(mumt); free(alphaphi); free(Bs2mualphaphi); free(Bs2mualphaphit);

}

      void updategam1(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p, double *newgam){ 
      int i,j;
      double maximum;
      double sumw = 0.0; 
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double *x= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *v1= malloc(sizeof(double)*N*p);
      double *v1t= malloc(sizeof(double)*p*N);
      double *wvg1= malloc(sizeof(double)*N);
      double *s= malloc(sizeof(double)*p);
      double *center= malloc(sizeof(double)*p);
      double *wvg2= malloc(sizeof(double)*N);
      double *Bs= malloc(sizeof(double)*N*p);
      double *v= malloc(sizeof(double)*N*p);
      double *Bs2= malloc(sizeof(double)*N*p);
      double *xBs2= malloc(sizeof(double)*N*p);
      double *xBs2t= malloc(sizeof(double)*p*N);
      double *dummy= malloc(sizeof(double)*N*p);
      double *dummy2= malloc(sizeof(double)*p*p);
      double *dummy3= malloc(sizeof(double)*p*p);
      double *dummyt= malloc(sizeof(double)*p*N);
      double *wty= malloc(sizeof(double)*N*p);
      double *wty2= malloc(sizeof(double)*N*p);
      double *F0t= malloc(sizeof(double)*p*p);
      double *F1t= malloc(sizeof(double)*p*p);
      double *F1tt= malloc(sizeof(double)*p*p);
      double *F2= malloc(sizeof(double)*p*p);
      double *cov1= malloc(sizeof(double)*p*p);
      double *cov11= malloc(sizeof(double)*p*p);
      double *cov1gam0= malloc(sizeof(double)*p*p);
      double *cov11gam0= malloc(sizeof(double)*p*p);
      double *cov2gam0= malloc(sizeof(double)*p*p);
      double *cov22gam0= malloc(sizeof(double)*p*p);
      double *cov2= malloc(sizeof(double)*p*p);
      double *cov22= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *A1= malloc(sizeof(double)*p*p);
      double *e2= malloc(sizeof(double)*N);
      double *max= malloc(sizeof(double)*N);
      double *e1= malloc(sizeof(double)*N);
      double *row= malloc(sizeof(double)*p);
      double *sign= malloc(sizeof(double)*p);
      double *u= malloc(sizeof(double)*p*p);
      double *ut= malloc(sizeof(double)*p*p);
      double *vtt= malloc(sizeof(double)*p*p);
      double *vttt= malloc(sizeof(double)*p*p);
      double *phimt= malloc(sizeof(double)*N*p);
      double *Bs2t= malloc(sizeof(double)*p*N);
      double *Bs2tphimt= malloc(sizeof(double)*p*N);
      double *mumt= malloc(sizeof(double)*N*p);
      double *alphaphi= malloc(sizeof(double)*N*p);
      double *Bs2mualphaphi= malloc(sizeof(double)*p*N);


      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=w[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs[i*p+j]=invW[i*p+j]*(double)1.0/phi0[j];
      
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2[j+i*p]=abc[1+i*3];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            phimt[i*p+j]=1.0/phi0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2t[j*N+i]=Bs2[j+i*p];
         for(i=0;i<N*p;i++)
            Bs2tphimt[i]=Bs2t[i]*phimt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs2[i*p+j]=Bs2tphimt[i+j*N];

      for(i=0;i<N;i++){
         wvg1[i]=w[i]*vg[0+i*2];
         wvg2[i]=w[i]*vg[1+i*2];
      }
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
            wty[i*p+j]=wvg2[i]*y[i*p+j];
            wty2[i*p+j]=wvg1[i]*y[i*p+j];
         }
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            dummy[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            dummyt[j*N+i]=dummy[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,dummyt,&N,&beta,dummy2,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs2[i*p+j]=x[i*p+j]*Bs2[i*p+j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            xBs2t[j*N+i]=xBs2[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty2,&p,xBs2t,&N,&beta,dummy3,&p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F0t[i*p+j]=dummy2[i*p+j]+dummy3[i*p+j];
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
             row[j]=Bs[i*p+j];
         }
         max[i]=maxi(row,p);
         e2[i]=max[i]*w[i]*vg[1+i*2];
      }
      
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
             row[j]=Bs2[i*p+j];
         }
         max[i]=maxi(row,p);
      }
      maximum=maxi(max,N);
      for(i=0;i<N;i++)
         e1[i]=maximum*w[i]*vg[0+i*2];
      for(i=0;i<N;i++)
            center[j]=0.0;
      Covariance(N,p,y,e2,center,cov1);
      Covariance(N,p,y,e1,center,cov11);
      Covariance(N,p,y,e1,center,cov2);
      Covariance(N,p,y,e2,center,cov22);
      double e2sum=0.0;
      double e1sum=0.0;
      for(i=0;i<N;i++){
          e1sum +=e1[i];
          e2sum +=e2[i];
      }
      for(i=0;i<p;i++){
         for(j=0;j<p;j++){
            cov1[i*p+j]=cov1[i*p+j]*e2sum;
            cov2[i*p+j]=cov2[i*p+j]*e1sum;
         }
      }
      for(i=0;i<p;i++){
         for(j=0;j<p;j++){
            cov11[i*p+j]=cov11[i*p+j]*e1sum;
            cov22[i*p+j]=cov22[i*p+j]*e2sum;
         }
      }
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0,&p,cov1,&p,&beta,cov1gam0,&p);
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0,&p,cov11,&p,&beta,cov11gam0,&p);
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0,&p,cov2,&p,&beta,cov2gam0,&p);
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0,&p,cov22,&p,&beta,cov22gam0,&p);
      if(e2sum==0){
        for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1t[i*p+j]=cov11gam0[i*p+j];
      }else if(e1sum==0){
        for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1t[i*p+j]=cov22gam0[i*p+j];
      }else{
       for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1t[i*p+j] = cov1gam0[i*p+j] +  cov2gam0[i*p+j];
      }
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1tt[i*p+j]=F1t[i+j*p];
    
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]=Bs[i*p+j]*mu0[j]+alpha0[j]/phi0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            mumt[i*p+j]=mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             alphaphi[i*p+j]=alpha0[j]/phi0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            Bs2t[j*N+i]=Bs2[j+i*p];
      for(i=0;i<N*p;i++)
         Bs2mualphaphi[i]=Bs2t[i]*mumt[i]+alphaphi[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v1[i*p+j]=Bs2mualphaphi[i+j*N];
      for(i=0;i<p;i++){
         for(j=0;j<N;j++){
            vt[i*N+j]=v[i+j*p];
            v1t[i*N+j]=v1[i+j*p];
         }
      }
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,vt,&N,&beta,A0,&p);
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty2,&p,v1t,&N,&beta,A1,&p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F2[i*p+j]=(F0t[i*p+j]-F1tt[i*p+j]-A0[i*p+j]-A1[i*p+j])/sumw;
      svd(p,p,F2,s,vtt,u);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            ut[i*p+j]=u[i+j*p];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            vttt[i*p+j]=vtt[i+j*p];
     dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,ut,&p,vttt,&p,&beta,newgam,&p);
      if(objgam(newgam,y,phi0,alpha0,mu0,w,vg,invW,abc,N,p) <objgam(gam0,y,phi0,alpha0,mu0,w,vg,invW,abc,N,p)){
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
             newgam[i*p+j]=gam0[i*p+j];
       }
       free(x); free(vt); free(v1); free(v1t); free(wvg1); free(s); 
       free(center); free(wvg2); free(Bs); free(v); free(Bs2); free(xBs2);
       free(xBs2t); free(dummy); free(dummy2); free(dummy3); free(dummyt);
       free(wty); free(wty2); free(F0t); free(F1t); free(F1tt); free(F2);
       free(cov1); free(cov1gam0); free(cov2gam0); free(cov2); free(A0);
       free(A1); free(e2); free(max); free(e1); free(row); free(sign);
       free(u); free(ut); free(vtt); free(phimt); free(Bs2t);
       free(Bs2tphimt); free(mumt); free(alphaphi); free(Bs2mualphaphi);
}








      double objgam(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p){

      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double *F0t= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *x= malloc(sizeof(double)*N*p);
      double *v= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *Bs= malloc(sizeof(double)*N*p);
      double *Bsmu= malloc(sizeof(double)*N*p);
      double *xBs= malloc(sizeof(double)*N*p);
      double *xBst= malloc(sizeof(double)*p*N);
      double *wtx= malloc(sizeof(double)*N*p);
      double *invUmat= malloc(sizeof(double)*N*p);
      double *wvg1= malloc(sizeof(double)*N);
      double *wvg2= malloc(sizeof(double)*N);
      
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
            invUmat[j+i*p]=abc[1+i*3];
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      double sumw=0;
      for(i=0;i<N;i++)
         sumw +=w[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs[i*p+j]=invW[i*p+j]*(double)1.0/phi0[j];
      for(i=0;i<N;i++)
         wvg2[i]=w[i]*vg[1+i*2];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
            wtx[j+i*p]=x[j+i*p]*wvg2[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             xBst[j*N+i]=xBs[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,xBst,&N,&beta,F0t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bsmu[i*p+j]=Bs[i*p+j]*mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]=Bsmu[i*p+j]+alpha0[j]/phi0[j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             vt[j*N+i]=v[j+i*p];
     dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,vt,&N,&beta,A0,&p); 
     double valMS= -1.0*(tr(F0t,p) -2.0*tr(A0,p))/sumw;

      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs[i*p+j]=invUmat[i*p+j]*(double)1.0/phi0[j];
      for(i=0;i<N;i++)
         wvg1[i]=w[i]*vg[0+i*2];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
            wtx[j+i*p]=x[j+i*p]*wvg1[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             xBst[j*N+i]=xBs[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,xBst,&N,&beta,F0t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bsmu[i*p+j]=Bs[i*p+j]*mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]=Bsmu[i*p+j]+alpha0[j]/phi0[j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             vt[j*N+i]=v[j+i*p];
     dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,vt,&N,&beta,A0,&p); 
     double valMUL= -1.0*(tr(F0t,p) -2.0*tr(A0,p))/sumw;
     double val = valMS + valMUL;

     free(F0t); free(A0); free(x); free(v); free(vt); free(Bs);
     free(Bsmu); free(xBs); free(xBst); free(wtx); free(invUmat);
     free(wvg1); free(wvg2);
     return val;
     
}


      double llik(double *x, double *mu, double *alpha, double *phi, double *cpl0, double **cpl, double **gam, double *wg, int N, int p, int G, double *pi, int *label){
      
      int g,j,i;
      double val=0.0;
      double cpl00[2];
      double wg0[2];
      double *logz= malloc(sizeof(double)*N*G);
      double *lablik= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *phi0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *zlog0= malloc(sizeof(double)*N);
      double *sum= malloc(sizeof(double)*N);
      int n=2;
      for(g=0;g<G;g++){
          for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             phi0[j]=phi[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
          for(i=0;i<n;i++){
              cpl00[i]= cpl0[g*n+i];
              wg0[i]= wg[g*n+i];
          }
          ddmsghyp(x,mu0,alpha0,phi0,cpl00,cpl[g],gam[g],wg0,zlog0,N,G,p);
          for(i=0;i<N;i++){
             logz[g+i*G]=zlog0[i];
           }  
      }
      for(i=0;i<N;i++)
         for(g=0;g<G;g++)
            lablik[i*G+g] =1.0;
      if(label)
         combinewk1(lablik,N,G,label);
      for(i=0;i<N;i++){
         sum[i]=0.0;
         for(g=0;g<G;g++){
            sum[i] +=(logz[i*G+g])*lablik[i*G+g]*pi[g];
         }
       }
      for(i=0;i<N;i++)
         val += log(sum[i]);
      free(logz); free(mu0); free(phi0); free(alpha0); free(zlog0); free(sum);
      free(lablik);
      return(val);
}

/*      void MAP(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double **sigma, double *wg, double *phi, int N, int G, int p, int v, double *pi, int *label, int *MAP){ 
      
      int i,g;
      double *w= malloc(sizeof(double)*N*G);
      double *C= malloc(sizeof(double)*G);
      v=1;
      weights(x,mu,alpha,cpl,cpl0,gam,sigma,wg,phi,N,G,p,v,pi,w);
      if(label)
        combinewk(w,N,G,label);
      for(i=0;i<N;i++){
         for(g=0;g<G;g++){
            C[g]=w[i*G+g];
         }
         MAP[i]=maxi_loc(C,G);
     free(w); free(C);
make sure you do not need to add 1 to MAP[i]
      }
}*/

     double tr(double *A, int N){
     int i;
     double trace=0.0;
     for(i=0;i<N;i++)
        trace +=A[i*N+i];
      return trace;
}

int maxi_loc(double *array, int size){
    int c;
    double max;
    max = array[0];
    int location = 0;
    for(c=0; c<size; c++){
        if(array[c] > max){
            max  = array[c];
            location = c;
        }
    }
    return location;
}

      void gig2(double *x, int N, int p, double *mu0, double *alpha0, double *cpl, double *cpl0, double *abc, double *phi0){
      int i,j;
      double sum=0.0;
      double v1,omega,a1;
      double *r= malloc(sizeof(double)*p);
      double *invS= malloc(sizeof(double)*p*p);
      double *b1= malloc(sizeof(double)*N);
      double *delta0= malloc(sizeof(double)*N);
      double *xmu= malloc(sizeof(double)*N*p);
      double *xmuphi= malloc(sizeof(double)*N);
      
         omega = cpl0[0];
      for(i=0; i<p; i++)
          sum += pow(alpha0[i],2)/phi0[i];
         a1=omega+sum;
      for(i=0; i<N; i++)
         for(j=0; j<p; j++)
            xmu[i*p+j]=x[i*p+j]-mu0[j];
      for(i=0; i<N; i++){
          xmuphi[i]=0.0;
          for(j=0; j<p; j++){
              xmuphi[i] += pow(xmu[i*p+j],2)*1/phi0[j];
          }
      }
      for(i=0; i<N; i++)
      b1[i] = omega + xmuphi[i]; 
      v1 = cpl0[1]-(double)p/(double)2.0;
      gig(b1,a1,v1,abc,N);
      free(r); free(invS); free(b1); free(delta0); free(xmu); free(xmuphi);
}

      void gig20(double *sr2sum, int N, int p, double *mu0, double *alpha0, double *cpl, double *cpl0, double *abc, double *phi0){
      int i;
      double sum=0.0;
      double v1,omega,a1;
      double *r= malloc(sizeof(double)*p);
      double *invS= malloc(sizeof(double)*p*p);
      double *b1= malloc(sizeof(double)*N);
      double *delta0= malloc(sizeof(double)*N);
      double *xmu= malloc(sizeof(double)*N*p);
      double *xmuphi= malloc(sizeof(double)*N);

         omega = cpl0[0];
      for(i=0; i<p; i++)
          sum += pow(alpha0[i],2)/phi0[i];
         a1=omega+sum;
//      for(i=0; i<N; i++)
//         for(j=0; j<p; j++)
//            xmu[i*p+j]=x[i*p+j]-mu0[j];
//      for(i=0; i<N; i++){
//          xmuphi[i]=0.0;
//          for(j=0; j<p; j++){
//              xmuphi[i] += pow(xmu[i*p+j],2)*1/phi0[j];
//          }
//      }
         
//      for(i=0; i<N; i++)
//      b1[i] = omega + xmuphi[i];
      for(i=0; i<N; i++)
          b1[i]=omega+sr2sum[i];
      v1 = cpl0[1]-(double)p/(double)2.0;
      gig(b1,a1,v1,abc,N);
      free(r); free(invS); free(b1); free(delta0); free(xmu); free(xmuphi);
}


      void gig(double *b, double a, double v, double *abc, int N){
      int i;

      double *sab= malloc(sizeof(double)*N);
      double *sba= malloc(sizeof(double)*N);
      double *kv1= malloc(sizeof(double)*N);
      double *kv= malloc(sizeof(double)*N);
      double *kv12= malloc(sizeof(double)*N);
      double *x= malloc(sizeof(double)*N);
      double *w= malloc(sizeof(double)*N);
      double *invw= malloc(sizeof(double)*N);
      double *logw= malloc(sizeof(double)*N);
     double *vec= malloc(sizeof(double)*N);

      for(i=0;i<N;i++)
          sab[i]=sqrt(a*b[i]);
      int exponscaled = 1;
      int nz = 1;
      double *r0= malloc(sizeof(double)*nz);
      int nSeq =1;
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
          BesselK(sab[i],nz,v+1,exponscaled,nSeq,r0,r1);
          kv1[i]=r0[0];
      }
      for(i=0;i<N;i++){
         BesselK(sab[i],nz,v,exponscaled,nSeq,r0,r1);
         kv[i]=r0[0];
      }
      for(i=0;i<N;i++)
         kv12[i]=kv1[i]/kv[i];
      for(i=0;i<N;i++)
          sba[i]=sqrt(b[i]/a);
      for(i=0;i<N;i++)
         w[i]=kv12[i]*sba[i];
      for(i=0;i<N;i++)
         invw[i]=kv12[i]/sba[i] -(double)2.0*v/b[i];
      for(i=0;i<N;i++)
         x[i]=v;
      double eps =1e-8;
      double d = 0.0001;
      int r = 6;
      int vv=2;
      grad2(x,N,sab,eps,d,r,vv,vec); 
      for(i=0;i<N;i++)
         logw[i]=log(sba[i])+vec[i];
      for(i=0;i<N;i++){
         abc[0+i*3]=w[i];
         abc[1+i*3]=invw[i];
         abc[2+i*3]=logw[i];
       }
      free(sab); free(sba); free(kv1); free(kv); free(kv12); free(x);
      free(w); free(invw); free(logw); free(vec); free(r1);
} 


      void mx_vec_mult(int n, int q, double *a, double *b, double *r){
      int i,k;
      for(i=0;i<n;i++){
         r[i]=0.0;
         for(k=0;k<q;k++){
            r[i] +=b[i*q+k]*a[k];
          }
       }
}

      void vec_mx_mult(int n, int q, double *a, double *b, double *r){
      int i,k;
      for(k=0;k<q;k++){
         r[k]=0.0;
         for(i=0;i<n;i++){
            r[k] +=a[i]*b[k+i*q];
         }
      }
} 




double maxi(double *array, int size){
    int c;
    double max;
    max = array[0];
    for(c=0; c<size; c++){
        if(array[c] > max){
            max  = array[c];
        }
    }
    return max;
}                             

      void ddghypFA(double *x, double *mu0, double *alpha0, double *cpl, double *sigma, double *val1, int N, int p,int log1){
     int i,j;
     int m=3;
     double omega, lambda, det;
     double pa;
     char notrans = 'N';
     double gamma = 1.0f;
     double beta = 0.0f;
     double sum=0;
     double dummy2=0.0;
     double lv[5];
     double *xmu  = malloc(sizeof(double)*N*p);
     double *dummy3  = malloc(sizeof(double)*N*p);
     double *dummy1  = malloc(sizeof(double)*p);
     double *dummy4  = malloc(sizeof(double)*N);
     double *lvx  = malloc(sizeof(double)*N*m);
     double *delta0= malloc(sizeof(double)*N);
     double *mx= malloc(sizeof(double)*N);
     double *kx= malloc(sizeof(double)*N);
     double *invS = malloc(sizeof(double)*p*p);
     omega = cpl[0];
     lambda= cpl[1];
     ginv(p,p,sigma,invS);
     vec_mx_mult(p,p,alpha0,invS,dummy1);
      for(i=0;i<p;i++)
      dummy2  +=dummy1[i]*alpha0[i];
      pa = omega + dummy2;
      mahalanobis(N,p,x,mu0,invS,delta0);
      for(i=0; i<N; i++){
          mx[i]= omega+delta0[i];
          kx[i]=sqrt(fabs(mx[i]*pa));
      }
      for(i=0; i<N; i++)
         for(j=0; j<p; j++)
            xmu[i*p+j] = x[i*p+j]-mu0[j];
      dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,invS,&p,xmu,&p,&beta,dummy3,&p);
      mx_vec_mult(N,p,alpha0,dummy3,dummy4);
      double nu=lambda -(double)p/(double)2.0;
      int exponscaled = 1;
      int nz = 1;
      int nSeq =1;
      double *r= malloc(sizeof(double)*N);
      double *r0= malloc(sizeof(double)*nz);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
        BesselK(kx[i],nz,nu,exponscaled,nSeq,r0,r1);
        r[i]=r0[0];
      }
      for(i=0; i<N; i++){
         lvx[0+i*m]=(lambda-(double)p/(double)2.0)*log(kx[i]);
         lvx[1+i*m]= log(r[i])-kx[i];
         lvx[2+i*m]= dummy4[i];
      }
     free(r); free(r1);
     nz =1;
     nu=lambda;
     double *r2= malloc(sizeof(double)*nz);
     double *r3= malloc(sizeof(double)*nz*nSeq);
     BesselK(omega,nz,nu,exponscaled,nSeq,r2,r3);
     determinant(sigma, p, p, &det);
     lv[0]=-((double)1.0/(double)2.0)*log(fabs(det))-((double)p/(double)2.0)*(log(2.0)+log(M_PI));

     lv[1]=omega-log(r2[0]);
     lv[2]=-(lambda/(double)2.0)*log(1);
     lv[3]=lambda*log(omega)*0;
     lv[4]=((double)p/(double)2-lambda)*log(fabs(pa));
     for(i=0;i<5;i++){
        sum = sum+lv[i];
     }
     for(i=0;i<N;i++){
        val1[i]=0;
        for(j=0;j<m;j++){
           val1[i] +=lvx[i*m+j];
        }
     val1[i]=val1[i]+sum;
     }
     free(xmu); free(dummy3); free(dummy1); free(dummy4); free(lvx);
     free(delta0); free(mx); free(kx); free(invS); free(r2); free(r3);

}

      void weightsFA(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p,int v, double *pi, double *w){  

      int g,i,j;
      int log1=1;
      double cpl1[2];
      double *zlog= malloc(sizeof(double)*N*G);
      double *mu0 = malloc(sizeof(double)*p);
      double *alpha0 = malloc(sizeof(double)*p);
      double *zlog0 = malloc(sizeof(double)*N);
      double sum =0.0;
      int n =2;
      for(g=0; g<G; g++){
          for(j=0; j<p; j++){
              mu0[j]=mu[g*p+j];
              alpha0[j] = alpha[g*p+j];
          }
          for(j=0;j<n;j++){
             cpl1[j]= cpl[g*n+j];
          }
   
          ddghypFA(x,mu0,alpha0,cpl1,sigma[g],zlog0,N,p,log1);
         for(j=0;j<N;j++){
             zlog[g+j*G] = zlog0[j];
         }
     }
     if(G>1){
        for(i=0;i<N;i++){
            sum=0.0;
            for(g=0;g<G;g++){
               w[i*G+g]=exp(v*(zlog[i*G+g]+log(pi[g])));
               if(w[i*G+g]==INFINITY) w[i*G+g]=1.0;
               sum +=w[i*G+g];
             }
    
         if(sum==0){
           for(g=0;g<G;g++)
              w[i*G+g]=(double)1.0/G;
         }else{
            for(g=0;g<G;g++){
               w[i*G+g]=w[i*G+g]/sum;
            }
      }
  }
  } else{
         for(i=0;i<N;i++)
            for(g=0;g<G;g++)
               w[i*G+g]=(double)1.0;
          }
}
   
  
      void gig2FA(double *x, int N, int p, double *mu0, double *alpha0, double *invS, double *cpl1, double *abc){
      int i;
      double sum=0.0;
      double omega, a1;
      double *r= malloc(sizeof(double)*p);
      double *delta0= malloc(sizeof(double)*N);
      double *b1= malloc(sizeof(double)*N);
      omega = exp(log(cpl1[0]));
      vec_mx_mult(p,p,alpha0,invS,r);
      for(i=0; i<p; i++)
          sum+=r[i]*alpha0[i];
      a1 = omega + sum;
      mahalanobis(N,p,x,mu0,invS,delta0);
      for(i=0;i<N;i++)
      for(i=0;i<N;i++)
         b1[i]=omega+delta0[i];
      double v1 = cpl1[1]-(double)p/(double)2.0;
      gigFA(b1,a1,v1,abc,N);
      free(r); free(delta0); free(b1);
      
}

      void gigFA(double *b, double a, double v, double *abc, int N){
      int i;

      double *sab= malloc(sizeof(double)*N);
      double *sba= malloc(sizeof(double)*N);
      double *kv1= malloc(sizeof(double)*N);
      double *kv= malloc(sizeof(double)*N);
      double *kv12= malloc(sizeof(double)*N);
      double *x= malloc(sizeof(double)*N);
      double *w= malloc(sizeof(double)*N);
      double *invw= malloc(sizeof(double)*N);
      double *logw= malloc(sizeof(double)*N);
      double *vec= malloc(sizeof(double)*N);

      for(i=0;i<N;i++)
          sab[i]=sqrt(fabs(a*b[i]));
      int exponscaled = 1;
      int nz = 1;
      double *r0= malloc(sizeof(double)*nz);
      int nSeq =1;
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
         BesselK(sab[i],nz,v+1,exponscaled,nSeq,r0,r1);
         kv1[i]=r0[0];
      }
      for(i=0;i<N;i++){
         BesselK(sab[i],nz,v,exponscaled,nSeq,r0,r1);
         kv[i]=r0[0];
      }
      for(i=0;i<N;i++)
         kv12[i]=kv1[i]/kv[i];
      for(i=0;i<N;i++)
          sba[i]=sqrt(fabs(b[i]/a));
      for(i=0;i<N;i++)
         w[i]=kv12[i]*sba[i];
      for(i=0;i<N;i++)
         invw[i]=kv12[i]/sba[i] -(double)2.0*v/b[i];
      for(i=0;i<N;i++)
         x[i]=v;
      double eps =1e-8;
      double d = 0.0001;
      int r = 6;
      int vv=2;
       
      grad(x,N,sab,eps,d,r,vv,vec); //what's vector come back
      for(i=0;i<N;i++)
         logw[i]=log(sba[i])-vec[i];
      for(i=0;i<N;i++){
         abc[0+i*3]=w[i];
         abc[1+i*3]=invw[i];
         abc[2+i*3]=logw[i];
      }
      free(sab); free(sba); free(kv1); free(kv); free(kv12); free(x);
      free(w); free(invw); free(logw); free(vec); free(r1);
}
     void EMgrstepFA(double *x, int N, int p, int G, int v, int q, int *label, double *mu, double *alpha, double *cpl, double **sigma, double **Lambda, double **err, double *pi){
     int i,j,g;
     int n=2;
     double cpl1[2];
     double cpl1new[2];
     double *wt= malloc(sizeof(double)*N);
     double *w= malloc(sizeof(double)*N*G);
     double *mu0= malloc(sizeof(double)*p);
     double *mu0new= malloc(sizeof(double)*p);
     double *alpha0= malloc(sizeof(double)*p);
     double *alpha0new= malloc(sizeof(double)*p);
        weightsFA(x,mu,alpha,cpl,sigma,N,G,p,v,pi,w);
      if(label)
        combinewk1(w,N,G,label);
      for(g=0;g<G;g++){
         for(j=0;j<p;j++){
            mu0[j]=mu[g*p+j];
            alpha0[j]=alpha[g*p+j];
         }
         for(j=0;j<n;j++){
            cpl1[j]=cpl[g*n+j];
          }
        for(i=0;i<N;i++){
           wt[i]=w[g+i*G];
         }
        updatemaScplM1(x,N,p,q,mu0,alpha0,mu0new,alpha0new,wt,cpl1,cpl1new,v,sigma[g]);
         for(j=0;j<p;j++){
            mu[g*p+j]=mu0new[j];
            alpha[g*p+j]=alpha0new[j];
         }
         for(j=0;j<n;j++){
             cpl[g*n+j]=cpl1new[j];
         }
        
        weightsFA(x,mu,alpha,cpl,sigma,N,G,p,v,pi,w);
        for(i=0;i<N;i++){
           wt[i]=w[g+i*G];
         }
        updatemaScplM2(x,N,p,q,mu0new,alpha0new,wt,cpl1new,sigma[g],v,Lambda[g],err[g]);
      }

      for(g=0;g<G;g++){
         pi[g]=0.0;
         for(i=0;i<N;i++){
            pi[g] +=w[g+i*G];
         }
         pi[g] /=N;
      }
      free(wt); free(w); free(mu0); free(mu0new); free(alpha0); free(alpha0new);
}




     void updatemaScplM1(double *x, int N, int p, int q, double *mu0, double *alpha0, double *mu0new, double *alpha0new, double *w, double *cpl1, double *cpl1new, int v, double *sigma ){

     int ncol=3;
     int i,j;
     double T =0.0;
     double sumw=0.0;
     double *abc= malloc(sizeof(double)*N*3);
     double *ABC= malloc(sizeof(double)*ncol);
     double *invS= malloc(sizeof(double)*p*p);
     double *u= malloc(sizeof(double)*N);
     double *t= malloc(sizeof(double)*N);
    
     if(!w)
       for(i=0;i<N;i++)
          w[i]=1;
     ginv(p,p,sigma,invS);
     gig2FA(x,N,p,mu0,alpha0,invS,cpl1,abc);
     for(i=0;i<N;i++)
         sumw +=w[i];
     for(j=0;j<ncol;j++){
        ABC[j]=0.0;
        for(i=0;i<N;i++){
            ABC[j] +=abc[j+i*ncol]*w[i];
        }
        ABC[j] /= sumw;
     }
     double  *alphaknown = NULL;
     if(!alphaknown){
        double A = ABC[0];
        double B = ABC[1];
        for(i=0;i<N;i++)
            u[i]=(B-abc[1+i*ncol])*w[i];
        for(i=0;i<N;i++){
            t[i]=(A*abc[1+i*ncol]-1)*w[i]; 
            T +=t[i];
        }
     for(j=0;j<p;j++){
        mu0new[j]=0.0;
        for(i=0;i<N;i++){
           mu0new[j] +=x[j+i*p]*t[i];
         }
         mu0new[j] /=T;
     }
     for(j=0;j<p;j++){
         alpha0new[j]=0.0;
        for(i=0;i<N;i++){
           alpha0new[j] +=x[j+i*p]*u[i];
        }
     alpha0new[j] /=T;
     }
   } else{

     for(j=0;j<p;j++)
        alpha0new[j]=alphaknown[j];
     for(j=0;j<p;j++){
        mu0new[j]=0.0;
        for(i=0;i<N;i++){
           mu0new[j] +=x[j+i*p]*abc[1+i*ncol]*w[i];
         }
           mu0new[j]=mu0new[j]-alpha0new[j]/ABC[1];
     }
  }
     for(j=0;j<p;j++)
        alpha0new[j]=v*alpha0new[j];
     double parol[2];
      
     parol[0]=cpl1[0];
     parol[1]=cpl1[1];
     updateolU(parol, ABC,2);
     updateolU(parol, ABC,2);
     cpl1new[0]=parol[0];
     cpl1new[1]=parol[1];
     free(abc); free(ABC); free(invS); free(u); free(t); 
}
     void updatemaScplM2(double *x, int N, int p, int q, double *mu0new, double *alpha0new, double *w,  double *cpl1new, double *sigma, int v, double *Lambda, double *err){

     char notrans = 'N';
     double alpha = 1.0f;
     double gama = 0.0f;
     int ncol=3;
     int i,j;
     double sumw=0.0;
     double *abc= malloc(sizeof(double)*N*3);
     double *ABC= malloc(sizeof(double)*ncol);
     double *A= malloc(sizeof(double)*p*p);
     double *r= malloc(sizeof(double)*p); 
     double *invS= malloc(sizeof(double)*p*p);
     double *wt= malloc(sizeof(double)*N);
     double *dummy= malloc(sizeof(double)*p*q);
     double *dummy4= malloc(sizeof(double)*q*q);
     double *dummy5= malloc(sizeof(double)*p*q);
     double *dummy6= malloc(sizeof(double)*p*p);
     double *inv= malloc(sizeof(double)*p*p);
     double *R= malloc(sizeof(double)*p*p);
     double *fi= malloc(sizeof(double)*p*p);
     double *fi1= malloc(sizeof(double)*p*p);
     double *fi1var= malloc(sizeof(double)*p*q);
     double *dia= malloc(sizeof(double)*q*q);
     double *var= malloc(sizeof(double)*p*p);
     double *vart= malloc(sizeof(double)*q*p);
     double *beta= malloc(sizeof(double)*q*p);
     double *vartfi1= malloc(sizeof(double)*q*p);
     double *vartfi1var= malloc(sizeof(double)*q*q);
     double *diavartfi1var= malloc(sizeof(double)*q*q);
     double *alphanewmat= malloc(sizeof(double)*N*p);
     double *alphanewmatt= malloc(sizeof(double)*p*N);
     double *abc1mat= malloc(sizeof(double)*N*p);
     double *abc2mat= malloc(sizeof(double)*N*p);
     double *x1= malloc(sizeof(double)*N*p);
     double *x2= malloc(sizeof(double)*N*p);
     double *x2t= malloc(sizeof(double)*p*N);
     double *term1= malloc(sizeof(double)*N*p);
     double *term1t= malloc(sizeof(double)*p*N);
     double *weightsmat= malloc(sizeof(double)*q*N);
     double *uhat= malloc(sizeof(double)*q*N);
     double *weightsmatuhat= malloc(sizeof(double)*N*q);
     double *weightsuhat= malloc(sizeof(double)*q*N);
     double *uhatsum= malloc(sizeof(double)*q);
     double *uhatt= malloc(sizeof(double)*N*q);
     double *uhatb= malloc(sizeof(double)*q*N);
     double *uhatbt= malloc(sizeof(double)*N*q);
     double *outer1= malloc(sizeof(double)*p*p);
     double *outer2= malloc(sizeof(double)*p*p);
     double *outer3= malloc(sizeof(double)*p*p);
     double *betavar= malloc(sizeof(double)*q*q);
     double *betaR= malloc(sizeof(double)*q*p);
     double *betaRbetat= malloc(sizeof(double)*q*q);
     double *betat= malloc(sizeof(double)*p*q);
     double *Iq= malloc(sizeof(double)*q*q);
     double *euu= malloc(sizeof(double)*q*q);
     double *euuinv= malloc(sizeof(double)*q*q);
     double *x2tuhatbt= malloc(sizeof(double)*p*q);
     double *alphanewmattuhatt= malloc(sizeof(double)*p*q);
     double *vanew= malloc(sizeof(double)*p*q);
     double *vanewt= malloc(sizeof(double)*q*p);
     double *fip1= malloc(sizeof(double)*p*p);
     double *alphauhat= malloc(sizeof(double)*p*q);
     double *alphauhatvart= malloc(sizeof(double)*p*p);
     double *vareuu= malloc(sizeof(double)*p*q);
     double *vareuuvart= malloc(sizeof(double)*p*p);
     double *fip2= malloc(sizeof(double)*p*p);
     double *vanewvanewt= malloc(sizeof(double)*p*p);
     double *finew= malloc(sizeof(double)*p);
     double *finewmat= malloc(sizeof(double)*p*p);

     if(!w)
       for(i=0;i<N;i++)
          w[i]=1;
     ginv(p,p,sigma,invS);
     gig2FA(x,N,p,mu0new,alpha0new,invS,cpl1new,abc);
     for(i=0;i<N;i++)
         sumw +=w[i];
     for(j=0;j<ncol;j++){
        ABC[j]=0.0;
        for(i=0;i<N;i++){
            ABC[j] +=abc[j+i*ncol]*w[i];
        }
        ABC[j] /= sumw;
     }
     for(i=0;i<p;i++)
         for(j=0;j<q;j++)
            var[i*q+j]=Lambda[i*q+j];
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
            fi[i*p+j]=err[i*p+j];
     
     double *fi2= malloc(sizeof(double)*p);
     double *fi3= malloc(sizeof(double)*p);
     for(i=0;i<p;i++)
         fi2[i]= err[i*p+i];
     double sumsq=0.0;
     for(i=0;i<p;i++)
        sumsq +=pow(fi2[i],2);
     for(i=0;i<p;i++)
        fi3[i]=fi2[i]/sumsq;
  
     for(i=0;i<p;i++){
        for(j=0;j<p;j++){
           fi1[i*p+j] = 0.0;
           if(i==j) fi1[i*p+i]=fi3[i];
        }
     }
     for(i=0;i<q;i++){
        for(j=0;j<q;j++){
           dia[i*q+j]=0.0;
           if(i==j) dia[i*q+i]=1.0;
        }
     }
     dgemm_(&notrans,&notrans,&q,&p,&p,&alpha,var,&q,fi1,&p,&gama,fi1var,&q);
     for(i=0;i<q;i++)
        for(j=0;j<p;j++)
           vart[i*p+j]=var[i+j*q];
     dgemm_(&notrans,&notrans,&p,&q,&p,&alpha,fi1,&p,vart,&p,&gama,vartfi1,&p);
     dgemm_(&notrans,&notrans,&q,&q,&p,&alpha,var,&q,vartfi1,&p,&gama,vartfi1var,&q);
     for(i=0;i<q;i++)
        for(j=0;j<q;j++)
            diavartfi1var[i*q+j]=dia[i*q+j]+vartfi1var[i*q+j];
     ginv(q,q,diavartfi1var,dummy4);
     dgemm_(&notrans,&notrans,&q,&p,&q,&alpha,dummy4,&q,fi1var,&q,&gama,dummy5,&q);
     dgemm_(&notrans,&notrans,&p,&p,&q,&alpha,vartfi1,&p,dummy5,&q,&gama,dummy6,&p);
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           inv[i*p+j]=fi1[i*p+j]-dummy6[i*p+j];
     dgemm_(&notrans,&notrans,&p,&q,&p,&alpha,inv,&p,vart,&p,&gama,beta,&p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            alphanewmat[i*p+j]=alpha0new[j];
        
     
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           abc1mat[j+i*p]=abc[0+i*ncol];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            x1[i*p+j]=x[i*p+j]-mu0new[j];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           term1[i*p+j]=x1[i*p+j]-alphanewmat[i*p+j]*abc1mat[i*p+j];
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           term1t[j*N+i] = term1[j+i*p];
     dgemm_(&notrans,&notrans,&N,&q,&p,&alpha,term1t,&N,beta,&p,&gama,uhat,&N);
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           abc2mat[j+i*p]=abc[1+i*ncol];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           term1[i*p+j]=x1[i*p+j]*abc2mat[i*p+j]-alphanewmat[i*p+j];
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           term1t[j*N+i] = term1[j+i*p];
     dgemm_(&notrans,&notrans,&N,&q,&p,&alpha,term1t,&N,beta,&p,&gama,uhatb,&N);
     for(i=0;i<N;i++)
        wt[i]=fabs(abc[1+i*ncol]*w[i]);
     Covariance(N,p,x,wt,mu0new,A);
   
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
            A[i*p+j]=A[i*p+j]*ABC[1];
     for(j=0;j<p;j++){
        r[j]=0.0;
        for(i=0;i<N;i++){
           r[j] +=x[j+i*p]*w[i];
        }
        r[j] /=sumw;
     }
           
     for(i=0;i<p;i++){
        for(j=0;j<p;j++){
           outer1[i*p+j]=alpha0new[j]*(r[i]-mu0new[i]);
           outer2[i*p+j]=(r[j]-mu0new[j])*alpha0new[i];
           outer3[i*p+j]=alpha0new[j]*alpha0new[i]*ABC[0];
        }
     }
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           R[i*p+j]=A[i*p+j]-(outer1[i*p+j]+outer2[i*p+j])+outer3[i*p+j];
     dgemm_(&notrans,&notrans,&q,&q,&p,&alpha,var,&q,beta,&p,&gama,betavar,&q);
     dgemm_(&notrans,&notrans,&p,&q,&p,&alpha,R,&p,beta,&p,&gama,betaR,&p);
     for(i=0;i<p;i++)
        for(j=0;j<q;j++)
           betat[i*q+j]=beta[i+j*p];
     dgemm_(&notrans,&notrans,&q,&q,&p,&alpha,betat,&q,betaR,&p,&gama,betaRbetat,&q);
     double sum =0.0;
     for(i=0;i<N;i++)
        sum +=abc[1+i*ncol]*w[i];
     for(i=0;i<q;i++){
        for(j=0;j<q;j++){
           Iq[i*q+j]=0.0;
           if(i==j)
             Iq[i*q+i]=sum;
        }
     }
     for(i=0;i<q;i++)
        for(j=0;j<q;j++)
           euu[i*q+j]=Iq[i*q+j]-sum*betavar[i*q+j]+betaRbetat[i*q+j]*sumw;
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           x2[i*p+j]=x[i*p+j]-mu0new[j];
     
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           x2[j+i*p]=x2[j+i*p]*w[i];
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           x2t[j*N+i]=x2[j+i*p];
     for(i=0;i<N;i++)
        for(j=0;j<q;j++)
           uhatbt[i*q+j]=uhatb[i+j*N];
     dgemm_(&notrans,&notrans,&q,&p,&N,&alpha,uhatbt,&q,x2t,&N,&gama,x2tuhatbt,&q);
    
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           alphanewmat[j+i*p]= alphanewmat[j+i*p]*w[i];
     for(j=0;j<p;j++)
        for(i=0;i<N;i++)
           alphanewmatt[j*N+i]=alphanewmat[j+i*p];
     for(i=0;i<N;i++)
        for(j=0;j<q;j++)
           uhatt[i*q+j]=uhat[i+j*N];
     dgemm_(&notrans,&notrans,&q,&p,&N,&alpha,uhatt,&q,alphanewmatt,&N,&gama,alphanewmattuhatt,&q);
     for(j=0;j<p;j++){
        for(i=0;i<q;i++){
            dummy[j*q+i]=x2tuhatbt[j*q+i]-alphanewmattuhatt[j*q+i];
        }
     }
     ginv(q,q,euu,euuinv);
     dgemm_(&notrans,&notrans,&q,&p,&q,&alpha,euuinv,&q,dummy,&q,&gama,vanew,&q);

     for(i=0;i<q;i++)
        for(j=0;j<p;j++)
           vart[i*p+j]=var[i+j*q];
     dgemm_(&notrans,&notrans,&p,&p,&q,&alpha,vart,&p,x2tuhatbt,&q,&gama,fip1,&p);
     for(i=0;i<q;i++)
        for(j=0;j<N;j++)
            weightsmat[i*N+j]=w[j];
     for(i=0;i<q*N;i++)
        weightsmatuhat[i]=weightsmat[i]*uhatt[i];
     for(i=0;i<q;i++)
        for(j=0;j<N;j++)
           weightsuhat[i*N+j]=weightsmatuhat[i+j*q];
     for(i=0;i<q;i++){
        uhatsum[i]=0.0;
        for(j=0;j<N;j++){
           uhatsum[i] += weightsuhat[i*N+j];
        }
     }
     for(i=0;i<p;i++)
        for(j=0;j<q;j++)
           alphauhat[i*q+j]=alpha0new[i]*uhatsum[j];
     dgemm_(&notrans,&notrans,&p,&p,&q,&alpha,vart,&p,alphauhat,&q,&gama,alphauhatvart,&p);
     dgemm_(&notrans,&notrans,&q,&p,&q,&alpha,euu,&q,var,&q,&gama,vareuu,&q);
     dgemm_(&notrans,&notrans,&p,&p,&q,&alpha,vart,&p,vareuu,&q,&gama,vareuuvart,&p);
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           fip2[i*p+j]=alphauhatvart[i*p+j]+vareuuvart[i*p+j];
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           R[i*p+j]=R[i*p+j]+(2.0/sumw)*(-fip1[i*p+j]+fip2[i*p+j]);
     for(i=0;i<p;i++)
        finew[i]=R[i*p+i];
     for(i=0;i<q;i++)
        for(j=0;j<p;j++)
           vanewt[i*p+j]=vanew[i+j*q];
     dgemm_(&notrans,&notrans,&p,&p,&q,&alpha,vanewt,&p,vanew,&q,&gama,vanewvanewt,&p);
     for(i=0;i<p;i++){
        for(j=0;j<p;j++){
           finewmat[i*p+j]=0;
           if(j==i)
              finewmat[i*p+i]=finew[i];
         }
     }
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
            sigma[i*p+j]=vanewvanewt[i*p+j]+finewmat[i*p+j];
     for(i=0;i<p;i++){
        for(j=0;j<q;j++){
           Lambda[i*q+j]=vanew[i*q+j];
        }
     }
     for(i=0;i<p;i++){
        for(j=0;j<p;j++){
           err[i*p+j]=0.0;
           if(i==j)
             err[i*p+i]=finew[i];
        }
    }
    free(abc); free(ABC); free(A); free(r); free(invS); free(wt); free(dummy);
    free(dummy4); free(dummy5); free(dummy6); free(inv); free(R); free(fi);
    free(fi1); free(fi1var); free(dia); free(var); free(vart); free(beta);
    free(vartfi1); free(vartfi1var); free(diavartfi1var); free(alphanewmat);
    free(alphanewmatt); free(abc1mat); free(abc2mat); free(x1); free(x2);
    free(x2t); free(term1); free(term1t); free(weightsmat); free(uhat); 
    free(weightsmatuhat); free(weightsuhat); free(uhatsum); free(uhatt);
    free(uhatb); free(uhatbt); free(outer1); free(outer2); free(outer3);
    free(betavar); free(betaR); free(betaRbetat); free(betat); free(Iq);
    free(euu); free(euuinv); free(x2tuhatbt); free(alphanewmattuhatt); free(vanew);
    free(vanewt); free(fip1); free(alphauhat); free(alphauhatvart); free(vareuu);
    free(vareuuvart); free(fip2); free(vanewvanewt); free(finew); free(finewmat); 

}


void MAPFA(double *x, double *mu, double *alpha, double *cpl, double **sigma,  int N, int G, int p, int v, double *pi, int *label, int *MAP){

      int i,g;
      double *w= malloc(sizeof(double)*N*G);
      double *C= malloc(sizeof(double)*G);
      v=1;
      weightsFA(x,mu,alpha,cpl,sigma,N,G,p,v,pi,w);
      if(label)
        combinewk1(w,N,G,label);
      for(i=0;i<N;i++){
         for(g=0;g<G;g++){
            C[g]=w[i*G+g];
         }
         MAP[i]=maxi_loc(C,G);
/*make sure you do not need to add 1 to MAP[i]*/
      }
      free(w); free(C);
}

      double llikFA(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int p, int G, double *pi){
      int n =2;
      int g,j,i;
      double cpl1[2];
      double val=0.0;
      double *logz= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *zlog0= malloc(sizeof(double)*N);
      double *sum= malloc(sizeof(double)*N);
      int log1=1;

      for(g=0;g<G;g++){
          for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
          for(j=0;j<n;j++){
            cpl1[j]=cpl[g*n+j];
          }
          ddghypFA(x,mu0,alpha0,cpl1,sigma[g], zlog0,N,p,log1);
          for(i=0;i<N;i++){
             logz[g+i*G]=zlog0[i];
           }
      }
      for(i=0;i<N;i++){
         sum[i]=0.0;
         for(g=0;g<G;g++){
            sum[i] +=exp(logz[i*G+g])*pi[g];
         }
       }
      for(i=0;i<N;i++)
         val += log(sum[i]);
      free(logz); free(mu0); free(alpha0); free(zlog0); free(sum);
      return(val);
}

/*void rgparGH(double *x, int N, int p, int G, double *w){
      int i,g;
      double *mu = malloc(sizeof(double)*G*p);
      double *mu0 = malloc(sizeof(double)*p);
      double *pi = malloc(sizeof(double)*G);
      double *wt = malloc(sizeof(double)*N);
      double *cpl = malloc(sizeof(double)*G*2);
      double *alpha = malloc(sizeof(double)*G*p);
      double **sigma    = malloc(sizeof(double*)*G);
      for(g=0; g<G; g++){
         sigma[g]      = malloc(sizeof(double)*p*p);
      }

      if(!w){
        for(i=0; i<N; i++)
           for(g=0; g<G; g++)
              w[i*G+g]= (double)1.0/(double)G;
      }
      for(g=0; g<G; g++){
        for(i=0; i<N; i++){
           wt[i]=w[g+i*G];
        }
      ipar(x,N,p,wt,G,mu0,alpha,sigma[g],cpl);
      for(i=0;i<p;i++){
         mu[g*p+i]=mu[i];
      }
}
      for(g=0; g<G; g++)
         pi[g]=(double)1.0/(double)G;
}*/

/*      void ipar(double *x, int N, int p, double *wt, int G, double *mu0, double *alpha0, double *sigma, double *cpl){
      int i, j,g;
      double *sampcov = malloc(sizeof(double)*p*p);
      double *alpha= malloc(sizeof(double)*G*p);
      double *vr= malloc(sizeof(double)*p*p);
      double *wr= malloc(sizeof(double)*p);
      double *var= malloc(sizeof(double)*p);
      double *wmean = malloc(sizeof(double)*p);

      if(!wt){    
        for(i=0; i<N; i++)
           wt[i] = (double)1.0;
      }    
      weightedmean(x,N,p,wt,wmean);
         for(j=0; j<p; j++)
         mu0[j] = randn(wmean[j],sqrt((double)1.0/(double)N));
      for(g=0; g<G; g++)
         for(j=0; j<p; j++)
            alpha[g*p+j] = randn(0.0,sqrt((double)1.0/(double)N));
      Covariance1(N, p, x, wt, mu0, sampcov);
      for(i=0;i<p;i++){
         for(j=0;j<p;j++){
            sigma[i*p+j]=0.0;
            if(i==j)
               sigma[i*p+i]=sampcov[i*p+i];
         }
      }
      for(i=0;i<p;i++){
         if(sigma[i*p+i]<0.1)
           sigma[i*p+i]=0.1;
      }
      eigen(p,sigma,wr,vr);
      for(i=0;i<p;i++){
         if(wr[i]<=0.0)
           variance(x,N,p,var);
           break;
      }
      for(j=0; j<p; j++)
          sigma[j*p+j]=var[j];
//    for(i=0;i<2;i++)
      for(g=0; g<G; g++){
          cpl[0+g*2]= 1.0;
          cpl[1+g*2]= -0.5;
      }

}*/
      void gig2GH(double *x, int N, int p, double *mu0, double *alpha0, double *invS, double *cpl1, double *abc){
      int i;
      double sum=0.0;
      double omega, a1;
      double *r= malloc(sizeof(double)*p);
      double *delta0= malloc(sizeof(double)*N);
      double *b1= malloc(sizeof(double)*N);
      omega = cpl1[0];
      vec_mx_mult(p,p,alpha0,invS,r);
      for(i=0; i<p; i++)
          sum+=r[i]*alpha0[i];
      a1 = omega + sum;
      mahalanobis(N,p,x,mu0,invS,delta0);
      for(i=0;i<N;i++)
         b1[i]=omega+delta0[i];
      double v1 = cpl1[1]-(double)p/(double)2.0;
      gigGH(b1,a1,v1,abc,N);
      free(r); free(delta0); free(b1);

}
      void gigGH(double *b, double a, double v, double *abc, int N){
      int i;

      double *sab= malloc(sizeof(double)*N);
      double *sba= malloc(sizeof(double)*N);
      double *kv1= malloc(sizeof(double)*N);
      double *kv= malloc(sizeof(double)*N);
      double *kv12= malloc(sizeof(double)*N);
      double *x= malloc(sizeof(double)*N);
      double *w= malloc(sizeof(double)*N);
      double *invw= malloc(sizeof(double)*N);
      double *logw= malloc(sizeof(double)*N);
      double *vec= malloc(sizeof(double)*N);

      for(i=0;i<N;i++)
          sab[i]=sqrt(a*b[i]);
      int exponscaled = 1;
      int nz = 1;
      int nSeq =1;
      double *r0= malloc(sizeof(double)*nz);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
         BesselK(sab[i],nz,v+1,exponscaled,nSeq,r0,r1);
         kv1[i]=r0[0];
      }
      for(i=0;i<N;i++){
        BesselK(sab[i],nz,v,exponscaled,nSeq,r0,r1);
        kv[i]=r0[0];
      }
      for(i=0;i<N;i++)
         kv12[i]=kv1[i]/kv[i];
      for(i=0;i<N;i++)
          sba[i]=sqrt(b[i]/a);
      for(i=0;i<N;i++)
         w[i]=kv12[i]*sba[i];
      for(i=0;i<N;i++)
         invw[i]=kv12[i]/sba[i] -(double)2.0*v/b[i];
      for(i=0;i<N;i++)
         x[i]=v;
      double eps =1e-8;
      double d = 0.0001;
      int r = 6;
      int vv=2;
      grad(x,N,sab,eps,d,r,vv,vec); 
      for(i=0;i<N;i++)
         logw[i]=log(sba[i])+vec[i];
//    for(g=0;g<3;g++)
      for(i=0;i<N;i++){
         abc[0+i*3]=w[i];
         abc[1+i*3]=invw[i];
         abc[2+i*3]=logw[i];
      }
     free(sab); free(sba); free(kv1); free(kv); free(kv12); free(x);
     free(w); free(invw); free(logw); free(vec); free(r1); free(r0);
}
     void updatemaScpl(double *x, int N, int p, double *mu0, double *alpha0,  double *w, double *cpl1, double *sigma, double v, double *mu0new, double *alpha0new, double *cpl1new){

     int ncol=3;
     int i,j;
     double T =0.0;
     double sumw=0.0;
     double parold[2];
     double parol[2];
     double *abc= malloc(sizeof(double)*N*3);
     double *ABC= malloc(sizeof(double)*ncol);
     double *A= malloc(sizeof(double)*p*p);
     double *invS= malloc(sizeof(double)*p*p);
     double *u= malloc(sizeof(double)*N);
     double *r= malloc(sizeof(double)*p);
     double *wt= malloc(sizeof(double)*N);
     double *t= malloc(sizeof(double)*N);
     double *R= malloc(sizeof(double)*p*p);
     double *dummy1= malloc(sizeof(double)*p*p);
     double *dummy2= malloc(sizeof(double)*p*p);
     double *dummy3= malloc(sizeof(double)*p*p);

     if(!w)
       for(i=0;i<N;i++)
          w[i]=1;
     ginv(p,p,sigma,invS);
     gig2GH(x,N,p,mu0,alpha0,invS,cpl1,abc);
     for(i=0;i<N;i++)
         sumw +=w[i];
     for(j=0;j<ncol;j++){
        ABC[j]=0.0;
        for(i=0;i<N;i++){
            ABC[j] +=abc[j+i*ncol]*w[i];
        }
        ABC[j] /=sumw;
     }
     double  *alphaknown = NULL;
     if(!alphaknown){
        double A = ABC[0];
        double B = ABC[1];
        for(i=0;i<N;i++)
            u[i]=(B-abc[1+i*ncol])*w[i];
        for(i=0;i<N;i++){
            t[i]=(A*abc[1+i*ncol]-1)*w[i];
            T +=t[i];
        }
     for(j=0;j<p;j++){
        mu0new[j]=0.0;
        for(i=0;i<N;i++){
           mu0new[j] +=x[j+i*p]*t[i];
         }
         mu0new[j] /=T;
     }
     for(j=0;j<p;j++){
         alpha0new[j]=0.0;
        for(i=0;i<N;i++){
           alpha0new[j] +=x[j+i*p]*u[i];
        }
     alpha0new[j] /=T;
     
     }
   } else{

     for(j=0;j<p;j++)
        alpha0new[j]=alphaknown[j];
     for(j=0;j<p;j++){
        mu0new[j]=0.0;
        for(i=0;i<N;i++){
           mu0new[j] +=x[j+i*p]*abc[1+i*ncol]*w[i];
         }
           mu0new[j]=mu0new[j]-alpha0new[j]/ABC[1];
     }
  }
     for(j=0;j<p;j++)
        alpha0new[j]=v*alpha0new[j];
      for(i=0;i<N;i++)
          wt[i]=abc[1+i*ncol]*w[i];
      Covariance(N,p,x,wt,mu0new,A);
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           A[i*p+j]=A[i*p+j]*ABC[1];
     for(j=0;j<p;j++){
        r[j]=0.0;
        for(i=0;i<N;i++){
            r[j] +=x[j+i*p]*w[i];
        }
        r[j] /=sumw;
      }
     for(i=0;i<p;i++){
        for(j=0;j<p;j++){
           dummy1[i*p+j]=(r[i]-mu0new[i])*alpha0new[j];
           dummy2[i*p+j]=alpha0new[i]*(r[j]-mu0new[j]);
           dummy3[i*p+j]=(alpha0new[i]*alpha0new[j])*ABC[0];
           R[i*p+j]=A[i*p+j]-(dummy1[i*p+j]+dummy2[i*p+j])+dummy3[i*p+j];
        }
     }
     for(i=0;i<p;i++){
        if(R[i*p+i]<0.00001)
          R[i*p+i]= 0.00001;
     }
     parold[0]=log(cpl1[0]);
     parold[1]=cpl1[1];
     parol[0]=exp(parold[0]);
     parol[1]=parold[1];
     int n=2;
     updateol(parol,ABC,n);
     updateol(parol,ABC,n);
     cpl1new[0]=parol[0];
     cpl1new[1]=parol[1];
     for(i=0;i<p;i++)
        for(j=0;j<p;j++)
           sigma[i*p+j]=R[i*p+j];
    free(abc); free(ABC); free(A); free(invS); free(u); free(r); free(wt);
    free(t); free(R); free(dummy1); free(dummy2); free(dummy3);
}
      void updateol(double *ol, double *ABC, int n){
      double Rp,Rn,f1,f2,bv;
      double lv1,lv0;
      int kmax =4;
      if(ABC[2] ==0){ 
         ol[1]=0;
      }else{
      double eps = 1e-8; 
      double d = 0.0001;
      int r =6;
      int vv=2;
      int N=1;
      double *vec=malloc(sizeof(double)*N);
      double *dummy1=malloc(sizeof(double)*N);
      double *dummy2=malloc(sizeof(double)*N);
      dummy1[0] = ol[1];
      dummy2[0] = ol[0];
       grad2(dummy1,N,dummy2,eps,d,r,vv,vec);
       bv = vec[0];
       ol[1]=ABC[2]*(ol[1]/bv);
       free(vec); free(dummy1); free(dummy2);
       } 
       int nz =1;
       double nu = ol[1];
//       int exponscaled = 1;
       int exponscaled = 0;
       int nSeq =1;
       double *r= malloc(sizeof(double)*nz);
       double *r2= malloc(sizeof(double)*nz);
       double *r1= malloc(sizeof(double)*nz*nSeq);
       BesselK(ol[0],nz,nu+1,exponscaled,nSeq,r,r1);
       BesselK(ol[0],nz,nu,exponscaled,nSeq,r2,r1);
       Rp = r[0]/r2[0];
       if(r[0]==INFINITY || r2[0]==INFINITY || r[0]==0 ||r2[0]==0){
          besselKnuAsym(ol[0],fabs(nu+1),kmax,&lv1);
          besselKnuAsym(ol[0],fabs(nu),kmax,&lv0);
       Rp = exp(lv1-lv0);
       }
 
       nu = -ol[1];
       BesselK(ol[0],nz,nu+1,exponscaled,nSeq,r,r1);
       BesselK(ol[0],nz,nu,exponscaled,nSeq,r2,r1);
       Rn=r[0]/r2[0];
       if(r[0]==INFINITY || r2[0]==INFINITY || r[0]==0 || r2[0]==0){
         besselKnuAsym(ol[0],fabs(nu+1),kmax,&lv1);
         besselKnuAsym(ol[0],fabs(nu),kmax,&lv0);
         Rn = exp(lv1-lv0);
       }
       f1 = Rp + Rn - (ABC[0]+ABC[1]);
       f2 = (pow(Rp,2)-(2*ol[1]+1.0)/ol[0]*Rp-1.0)+(pow(Rn,2)-(2.0*(-1*ol[1])+1.0)/ol[0]*Rn-1.0);
       if((ol[0]-f1/f2)>0){
         ol[0]=ol[0]-f1/f2;}
      free(r); free(r2); free(r1);
}
      void ddghypGH(double *x, double *mu0, double *alpha0, double *cpl, double *sigma, double *val1, int N, int p,int log1){
     int i,j;
     int m=3;
     double omega, lambda, det;
     double pa;
     char notrans = 'N';
     double gamma = 1.0f;
     double beta = 0.0f;
     double sum=0;
     double dummy2=0.0;
     double lv[5];
     double *xmu  = malloc(sizeof(double)*N*p);
     double *dummy3  = malloc(sizeof(double)*N*p);
     double *dummy1  = malloc(sizeof(double)*p);
     double *dummy4  = malloc(sizeof(double)*N);
     double *lvx  = malloc(sizeof(double)*N*m);
     double *delta0= malloc(sizeof(double)*N);
     double *mx= malloc(sizeof(double)*N);
     double *kx= malloc(sizeof(double)*N);
     double *invS = malloc(sizeof(double)*p*p);
     omega = cpl[0];
     lambda= cpl[1];
     ginv(p,p,sigma,invS);
     vec_mx_mult(p,p,alpha0,invS,dummy1);
      for(i=0;i<p;i++)
      dummy2  +=dummy1[i]*alpha0[i];
      pa = omega + dummy2;
      mahalanobis(N,p,x,mu0,invS,delta0);
      for(i=0; i<N; i++){
          mx[i]= omega+delta0[i];
          kx[i]=sqrt(mx[i]*pa);
      }
      for(i=0; i<N; i++)
         for(j=0; j<p; j++)
            xmu[i*p+j] = x[i*p+j]-mu0[j];
      dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,invS,&p,xmu,&p,&beta,dummy3,&p);
      mx_vec_mult(N,p,alpha0,dummy3,dummy4);
      double nu=lambda -(double)p/(double)2.0;
      int exponscaled = 1;
      int nz = 1;
      int nSeq =1;
      double *r= malloc(sizeof(double)*N);
      double *r0= malloc(sizeof(double)*nz);
      double *r1= malloc(sizeof(double)*nz*nSeq);
      for(i=0;i<N;i++){
         BesselK(kx[i],nz,nu,exponscaled,nSeq,r0,r1);
         r[i]=r0[0];
      }

      for(i=0; i<N; i++){
         lvx[0+i*m]=(lambda-(double)p/(double)2.0)*log(kx[i]);
         lvx[1+i*m]= log(r[i])-kx[i];
         lvx[2+i*m]= dummy4[i];
      }
     free(r); free(r0); free(r1);
     nz =1;
     nu=lambda;
     double *r2= malloc(sizeof(double)*nz);
     double *r3= malloc(sizeof(double)*nz*nSeq);
     BesselK(omega,nz,nu,exponscaled,nSeq,r2,r3);
     determinant(sigma, p, p, &det);
     lv[0]=-((double)1.0/(double)2.0)*log(det)-((double)p/(double)2.0)*(log(2.0)+log(M_PI));
     lv[1]=omega-log(r2[0]);
     lv[2]=-(lambda/(double)2.0)*log(1);
     lv[3]=lambda*log(omega)*0;
     lv[4]=((double)p/(double)2-lambda)*log(pa);
     for(i=0;i<5;i++)
        sum = sum+lv[i];
     for(i=0;i<N;i++){
        val1[i]=0;
        for(j=0;j<m;j++){
           val1[i] +=lvx[i*m+j];
        }
     val1[i]=val1[i]+sum;
     }
     if(!log1){
       for(i=0;i<N;i++)
          val1[i]=exp(val1[i]);
     }
     free(xmu); free(dummy3); free(dummy1); free(dummy4); free(lvx);
     free(delta0); free(mx); free(kx); free(invS); free(r2); free(r3);
}

      void weightsGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p,int v, double *pi, double *w){

      int g,i,j;
      int log1=1;
      int n =2;
      double cpl1[2];
      double *zlog= malloc(sizeof(double)*N*G);
      double *mu0 = malloc(sizeof(double)*p);
      double *alpha0 = malloc(sizeof(double)*p);
      double *C = malloc(sizeof(double)*G);
      double *max = malloc(sizeof(double)*N);
      double *zlog0 = malloc(sizeof(double)*N);
      double sum =0.0;
      for(g=0; g<G; g++){
          for(j=0; j<p; j++){
              mu0[j]=mu[g*p+j];
              alpha0[j] = alpha[g*p+j];
          }
         for(j=0;j<n;j++){
            cpl1[j]=cpl[g*n+j];
          }
          ddghypGH(x,mu0,alpha0,cpl1,sigma[g],zlog0,N,p,log1);
         for(j=0;j<N;j++){
             zlog[g+j*G] = zlog0[j];
         }
     }
     if(G>1){
        for(i=0;i<N;i++){
            for(g=0;g<G;g++){
               C[g]=v*(zlog[i*G+g]+log(pi[g]));
            }
            max[i]=maxi(C,G);
        }
        for(i=0;i<N;i++){
            sum=0.0;
            for(g=0;g<G;g++){
               w[i*G+g]=v*(zlog[i*G+g]+log(pi[g]))-max[i];
               w[i*G+g]=exp(w[i*G+g]);
               sum +=w[i*G+g];
             }
            for(g=0;g<G;g++){
               w[i*G+g]= w[i*G+g]/sum;
            }
         if(sum==0){
           for(g=0;g<G;g++)
              w[i*G+g]=(double)1.0/G;
         }
      }
   }else{
         for(i=0;i<N;i++)
            for(g=0;g<G;g++)
               w[i*G+g]=(double)1.0;
          }
     free(zlog); free(mu0); free(alpha0); free(C); free(max); free(zlog0);
}
      void EMgrstepGH (double *x, int N, int p, int G, int v, int *label, double *mu, double *alpha, double *cpl, double **sigma, double *pi ){
     int i,j,g;
     int n =2;
     double cpl1[2];
     double cpl1new[2];
     double *wt= malloc(sizeof(double)*N);
     double *w= malloc(sizeof(double)*N*G);
     double *mu0= malloc(sizeof(double)*p);
     double *mu0new= malloc(sizeof(double)*p);
     double *alpha0= malloc(sizeof(double)*p);
     double *alpha0new= malloc(sizeof(double)*p);
        weightsGH(x,mu,alpha,cpl,sigma,N,G,p,v,pi,w);
      if(label)
        combinewk1(w,N,G,label);
      for(g=0;g<G;g++){
         for(j=0;j<p;j++){
            mu0[j]=mu[g*p+j];
            alpha0[j]=alpha[g*p+j];
         }
         for(j=0;j<n;j++){
            cpl1[j]=cpl[g*n+j];
          }
        for(i=0;i<N;i++){
           wt[i]=w[g+i*G];
         }
        updatemaScpl(x,N,p,mu0,alpha0,wt,cpl1,sigma[g],v,mu0new,alpha0new,cpl1new);
         for(j=0;j<p;j++){
            mu[g*p+j]=mu0new[j];
            alpha[g*p+j]=alpha0new[j];
         }
         for(j=0;j<n;j++){
             cpl[g*n+j]=cpl1new[j];
         }
      }
          
      for(g=0;g<G;g++){
         pi[g]=0.0;
         for(i=0;i<N;i++){
            pi[g] +=w[g+i*G];
         }
         pi[g] /=N;
      }
      free(wt); free(mu0); free(mu0new); free(alpha0); free(alpha0new);
}


      void MAPGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p, int v, double *pi, int *label, int *MAP){
        
      int i,g;
      double *w= malloc(sizeof(double)*N*G);
      double *C= malloc(sizeof(double)*G);
      v=1;
      weightsGH(x,mu,alpha,cpl,sigma,N,G,p,v,pi,w);
      if(label)
        combinewk(w,N,G,label);
      for(i=0;i<N;i++){
         for(g=0;g<G;g++){
            C[g]=w[i*G+g];
         } 
         MAP[i]=maxi_loc(C,G);
// make sure you do not need to add 1 to MAP[i]
     }
     free(w); free(C);
}

      double llikGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int p, int G, double *pi){
      int n =2;
      int g,j,i;
      double cpl1[2];
      double val=0.0;
      double *logz= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *zlog0= malloc(sizeof(double)*N);
      double *sum= malloc(sizeof(double)*N);
      int log1=1;

      for(g=0;g<G;g++){
          for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
          for(j=0;j<n;j++){
            cpl1[j]=cpl[g*n+j];
          }
          ddghypGH(x,mu0,alpha0,cpl1,sigma[g], zlog0,N,p,log1);
          for(i=0;i<N;i++){
             logz[g+i*G]=zlog0[i];
           }
      }
      for(i=0;i<N;i++){
         sum[i]=0.0;
         for(g=0;g<G;g++){
            sum[i] +=exp(logz[i*G+g])*pi[g];
         }
       }
      for(i=0;i<N;i++)
         val += log(sum[i]);
      free(logz); free(mu0); free(alpha0); free(zlog0); free(sum);
      return(val);
}



      void dmsghypMS(double *y, double *mu0, double *alpha0, double *cpl, double *gam, double *sigma0, double *val2, int N, int p){
     int i,j;
     int nz;
     char notrans = 'N';
     double gamma = 1.0f;
     double beta = 0.0f;
     double *x    = malloc(sizeof(double)*N*p);
     double *m0    = malloc(sizeof(double)*N*p);
     double *mx    = malloc(sizeof(double)*N*p);
     double *kx    = malloc(sizeof(double)*N*p);
     double *lx1    = malloc(sizeof(double)*N*p);
     double *lx2    = malloc(sizeof(double)*N*p);
     double *matrix    = malloc(sizeof(double)*N*p);
     double *sum    = malloc(sizeof(double)*N);
     double *lvsum    = malloc(sizeof(double)*p);
     double *lx3    = malloc(sizeof(double)*N*p);
     double *lv1    = malloc(sizeof(double)*p);
     double *lv2   = malloc(sizeof(double)*p);
     double *lv3    = malloc(sizeof(double)*p);
     double *pa    = malloc(sizeof(double)*p);
     double *row    = malloc(sizeof(double)*p);
     double *nu1    = malloc(sizeof(double)*p);
     double *chi    = malloc(sizeof(double)*p);
     double *psi    = malloc(sizeof(double)*p);
     double *lambda    = malloc(sizeof(double)*p);
     double *dummy  = malloc(sizeof(double)*N*p);
     double *yy  = malloc(sizeof(double)*N*p);
     double *xmu  = malloc(sizeof(double)*N*p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           yy[i*p+j]= y[i+j*N];
     dgemm_(&notrans,&notrans,&p,&N,&p,&gamma,gam,&p,yy,&p,&beta,x,&p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            xmu[i*p+j] = x[i*p+j]-mu0[j];
            
     for(i=0;i<p;i++){
        chi[i]=cpl[0+i*2];
        psi[i]=cpl[0+i*2];
        lambda[i]=cpl[1+i*2];
     }
     for(i=0;i<p;i++)
        pa[i]=psi[i]+pow(alpha0[i],2.0)/sigma0[i];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            m0[i*p+j]=pow(xmu[i*p+j],2.0)*(double)1.0/sigma0[j];
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           mx[i*p+j]=m0[i*p+j]+chi[j];
     for(i=0;i<N;i++){
        for(j=0;j<p;j++){
           kx[i*p+j]=sqrt(mx[i*p+j]*pa[j]);
        }
     }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           dummy[i*p+j]=log(mx[i*p+j])-log(pa[j]);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           lx1[i*p+j]=dummy[i*p+j]*(lambda[j]-(double)1.0/(double)2.0)/(double)2.0;
     int exponscaled=1;
     nz = 1;
     int nSeq = 1;
     for(j=0;j<p;j++)
        nu1[j] = lambda[j] -(double)1/(double)2;
      double *r0= malloc(sizeof(double)*nz);
      double *r= malloc(sizeof(double)*p);
      double *r1= malloc(sizeof(double)*nz*nSeq);

     for(i=0;i<N;i++){
        for(j=0;j<p;j++){
           row[j]=kx[i*p+j];
           BesselK(row[j],nz,nu1[j],exponscaled,nSeq,r0,r1);
           r[j]=r0[0];
        }
        for(j=0;j<p;j++){
           lx2[i*p+j]=log(r[j])-kx[i*p+j];
        }
     }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
            lx3[i*p+j]=xmu[i*p+j]*alpha0[j]/sigma0[j];
     for(i=0;i<p;i++)
        lv1[i]=-(double).5*log(sigma0[i])-(double).5*(log((double)2)+log(M_PI));
     for(j=0;j<p;j++){
        BesselK(chi[j],nz,lambda[j],exponscaled,nSeq,r0,r1);
        r[j]=r0[0];
     }
     for(i=0;i<p;i++){
        lv2[i]=-(log(r[i])-chi[i]);
         lv3[i]=0;
     }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           matrix[i*p+j]=lx1[i*p+j]+lx2[i*p+j]+lx3[i*p+j];
     for(i=0;i<N;i++){
         sum[i]=0.0;
         for(j=0;j<p;j++){
            sum[i] +=matrix[i*p+j];
         }
      }
      double lvtot=0.0;
      for(j=0;j<p;j++){
         lvsum[j]=lv1[j]+lv2[j]+lv3[j];
         lvtot +=lvsum[j];
      }
     for(i=0;i<N;i++){
        val2[i]= sum[i]+lvtot;
    }

      free(x); free(m0); free(mx); free(kx); free(lx1); free(lx2); 
      free(matrix); free(sum); free(lvsum); free(lx3); free(lv1);
      free(lv2); free(lv3); free(pa); free(row); free(nu1); free(chi); 
      free(psi); free(lambda); free(dummy); free(xmu); free(r0);
      free(r); free(r1);
}

      void weightsMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *phi, int N, int G, int p,int v, double *pi, double *w){

      int g,i,j;
      double *zlog= malloc(sizeof(double)*N*G);
      double *mu0 = malloc(sizeof(double)*p);
      double *alpha0 = malloc(sizeof(double)*p);
      double *sigma0 = malloc(sizeof(double)*p);
      double *zlog0 = malloc(sizeof(double)*N);
      double sum =0.0;
      for(g=0; g<G; g++){
          for(j=0; j<p; j++){
              mu0[j]=mu[g*p+j];
              alpha0[j] = alpha[g*p+j];
              sigma0[j]= phi[g*p+j];
          }
          dmsghypMS(x,mu0,alpha0,cpl[g],gam[g],sigma0,zlog0,N,p);
         for(j=0;j<N;j++){
             zlog[g+j*G] = zlog0[j];
         }
     }
     if(G>1){
        for(i=0;i<N;i++){
            sum=0.0;
            for(g=0;g<G;g++){
               w[i*G+g]=exp(v*(zlog[i*G+g]+log(pi[g])));
               sum +=w[i*G+g];
             }
            for(g=0;g<G;g++){
               w[i*G+g] /=sum;
            }
         if(sum==0){
           for(g=0;g<G;g++)
              w[i*G+g]=(double)1.0/G;
         }
      }
   }else{
         for(i=0;i<N;i++){
            for(g=0;g<G;g++){
               w[i*G+g]=(double)1.0;
            }
          }
        }
     free(zlog); free(mu0); free(alpha0); free(sigma0); free(zlog0); 
}


       void EMgrstepMS(double *y, double **gam, double *mu, double *phi, double *alpha, double **cpl, int N, int p, int G, double *pi, int v, int *label,int m) {
       int i,j,g;
      double *w= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *phi0= malloc(sizeof(double)*p);
      double *mu0new= malloc(sizeof(double)*p);
      double *alpha0new= malloc(sizeof(double)*p);
      double *phi0new= malloc(sizeof(double)*p);
      double *cplnew= malloc(sizeof(double)*p*2);
      double *newgam= malloc(sizeof(double)*p*p);
      double *wt= malloc(sizeof(double)*N);
      if(label)
        combinewk1(w,N,G,label);
      for(g=0;g<G;g++){
         weightsMS(y,mu,alpha,cpl,gam,phi,N,G,p,v,pi,w);
         if(label)
           combinewk1(w,N,G,label);
         for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             phi0[j]=phi[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
      for(i=0;i<N;i++){
         wt[i]=w[g+i*G];
      }
      updatemaScplpMS(y,gam[g],mu0,phi0,alpha0,cpl[g],wt,N,p,mu0new,alpha0new,phi0new,newgam,cplnew,v,m);
      for(j=0;j<p;j++){
          mu[g*p+j]=mu0new[j];
          phi[g*p+j]=phi0new[j];
          alpha[g*p+j]=alpha0new[j];
      }
      for(j=0;j<p*2;j++)
          cpl[g][j]=cplnew[j];
      for(j=0;j<p*p;j++)
          gam[g][j]=newgam[j];
    }
      for(g=0;g<G;g++){
         pi[g]=0.0;
         for(i=0;i<N;i++){
             pi[g] +=w[g+i*G];
         }
         pi[g] /=N;
      }
      free(w); free(mu0); free(alpha0); free(phi0); free(mu0new);
      free(alpha0new); free(phi0new); free(cplnew); free(newgam); free(wt);

} 
       void updatemaScplpMS(double *y, double *gam, double *mu0, double *sigma0, double *alpha0 , double *cpl, double *wt, int N, int p, double *mu0new, double *alpha0new, double *sigma0new, double *newgam, double *cplnew, int v, int m){
      int i,j;
     char notrans = 'N';
     double alpha = 1.0f;
     double beta = 0.0f;
      double sumw = 0.0;
      double *x= malloc(sizeof(double)*N*p);
      double *sr2= malloc(sizeof(double)*N*p);
      double *xmu= malloc(sizeof(double)*N*p);
      double *W= malloc(sizeof(double)*N*p);
      double *invW= malloc(sizeof(double)*N*p);
      double *logW= malloc(sizeof(double)*N*p);
      double *u= malloc(sizeof(double)*N*p);
      double *t= malloc(sizeof(double)*N*p);
      double *A= malloc(sizeof(double)*p);
      double *B= malloc(sizeof(double)*p);
      double *C= malloc(sizeof(double)*p);
      double *T= malloc(sizeof(double)*p);
      double *xt= malloc(sizeof(double)*p);
      double *xu= malloc(sizeof(double)*p);
      double *Ax= malloc(sizeof(double)*p);
      double *ax= malloc(sizeof(double)*p);
      double *omega= malloc(sizeof(double)*p);
      double *test1= malloc(sizeof(double)*p*2);
      double *test2= malloc(sizeof(double)*p*3);
      double *row1= malloc(sizeof(double)*2);
      double *row2= malloc(sizeof(double)*3);
      double *yy= malloc(sizeof(double)*N*p);
      

      if(!wt){
        for(i=0;i<N;i++)
           wt[i]=1;
      }
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           yy[i*p+j]= y[i+j*N];
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam,&p,yy,&p,&beta,x,&p);
     for(i=0; i<N;i++){
         for(j=0;j<p;j++){
            x[i*p+j]=pow(x[i*p+j]-mu0[j],2);
            sr2[i*p+j]=x[i*p+j]*1/sigma0[j];
         }
      }
      gig2pMS(sr2,N,p,cpl,sigma0,alpha0,W,invW,logW);
      if(m%2==0) 
         updategam2MS(gam,yy,sigma0,alpha0,mu0,wt,invW,newgam,N,p);
      else
         updategam1MS(gam,yy,sigma0,alpha0,mu0,wt,invW,newgam,N,p);
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,newgam,&p,yy,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=wt[i];      
      for(j=0;j<p;j++){
         A[j]=0.0;
         for(i=0; i<N;i++){
             A[j] +=W[j+i*p]*wt[i];
         }
         A[j] /=sumw;
      }
      for(j=0;j<p;j++){
          B[j]=0.0;
         for(i=0; i<N;i++){
            B[j] +=invW[j+i*p]*wt[i];
         }  
         B[j] /=sumw;
      }
      for(j=0;j<p;j++){
          C[j]=0.0;
         for(i=0; i<N;i++){
            C[j] +=logW[j+i*p]*wt[i];
         }
         C[j] /=sumw;
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            u[i*p+j]= (B[j]-invW[i*p+j])*wt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            t[i*p+j]=(invW[i*p+j]*A[j]-1.0)*wt[i];
      for(j=0;j<p;j++){
         T[j]=0.0;
         for(i=0;i<N;i++){
            T[j] +=t[j+i*p];       
         }
      }
      for(j=0;j<p;j++){
         xt[j]=0.0;
         for(i=0;i<N;i++){
            xt[j] +=x[j+i*p]*t[j+i*p];
         }
      }
      for(j=0;j<p;j++)
         mu0new[j]=xt[j]/T[j];
      for(j=0;j<p;j++){
         xu[j]=0.0;
         for(i=0;i<N;i++){
            xu[j] +=x[j+i*p]*u[j+i*p];
         }
      }
      for(j=0;j<p;j++)
         alpha0new[j]=v*xu[j]/T[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             xmu[i*p+j]=pow(x[i*p+j]-mu0new[j],2)*invW[i*p+j];
      for(j=0;j<p;j++){
         Ax[j]=0.0;
         for(i=0;i<N;i++){
            Ax[j] +=xmu[j+i*p]*wt[i];
          }
          Ax[j] /=sumw;
       }
      for(j=0;j<p;j++){
         ax[j]=0.0;
         for(i=0;i<N;i++){
            ax[j] +=x[j+i*p]*wt[i];
         }
         ax[j] /=sumw;
      }
      for(j=0;j<p;j++)
         sigma0new[j]=Ax[j]-2.0*(ax[j]-mu0new[j])*alpha0new[j]+pow(alpha0new[j],2)*A[j];
      for(j=0;j<p;j++)
         omega[j]=cpl[0+j*2];

/*   for(i=0;i<2;i++)*/
       for(j=0;j<p;j++){
           test1[0+j*2]=omega[j];
           test1[1+j*2]=cpl[1+j*2];
       }
       for(j=0;j<p;j++){
          test2[0+j*3]=A[j];
          test2[1+j*3]=B[j];
          test2[2+j*3]=C[j];
       }
       for(j=0;j<p;j++){
          for(i=0;i<2;i++){
             row1[i]=test1[j*2+i];
          }
          for(i=0;i<3;i++){
             row2[i]=test2[j*3+i];
          }
          updateol(row1,row2,2);
          updateol(row1,row2,2);
          for(i=0;i<2;i++){
             cplnew[j*2+i]=row1[i];
          }
       }
       free(x); free(sr2); free(xmu); free(W); free(invW); free(logW);
       free(u); free(t); free(A); free(B); free(C); free(T); free(xt);
       free(xu); free(Ax); free(ax); free(omega); free(test1); free(test2);
       free(row1); free(row2);
}

      double llikMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *sigma, int N, int p, int G, double *pi){
      int g,j,i; 
      double val=0.0;
      double *logz= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *sigma0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *zlog0= malloc(sizeof(double)*N);
      double *sum= malloc(sizeof(double)*N);
           
      for(g=0;g<G;g++){
          for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             sigma0[j]=sigma[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
          dmsghypMS(x,mu0,alpha0,cpl[g],gam[g],sigma0,zlog0,N,p);
          for(i=0;i<N;i++){
             logz[g+i*G]=zlog0[i];
           }
      }
      for(i=0;i<N;i++){
         sum[i]=0.0;
         for(g=0;g<G;g++){
            sum[i] +=exp(logz[i*G+g])*pi[g];
         }
       }
      for(i=0;i<N;i++)
         val += log(sum[i]);
      free(logz); free(mu0); free(sigma0); free(alpha0); free(zlog0); free(sum);
      return(val);
}
      
      void MAPMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *sigma, int N, int G, int p, int v, double *pi, int *label, int *MAP){

      int i,g;
      double *w= malloc(sizeof(double)*N*G);
      double *C= malloc(sizeof(double)*G);
      v=1;
      weightsMS(x,mu,alpha,cpl,gam,sigma,N,G,p,v,pi,w);
      if(label)
        combinewk(w,N,G,label);
      for(i=0;i<N;i++){
         for(g=0;g<G;g++){
            C[g]=w[i*G+g];
         }
         MAP[i]=maxi_loc(C,G);
/*make sure you do not need to add 1 to MAP[i]*/
      }
      free(w); free(C);
}

      void updategam2MS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *wt, double *invW, double *newgam, int N, int p){
      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double sumw = 0.0;
      double *x= malloc(sizeof(double)*N*p);
      double *yy= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *v= malloc(sizeof(double)*N*p);
      double *Bs= malloc(sizeof(double)*N*p);
      double *wty= malloc(sizeof(double)*N*p);
      double *xBs= malloc(sizeof(double)*N*p);
      double *Bsmu= malloc(sizeof(double)*N*p);
      double *xBst= malloc(sizeof(double)*p*N);
      double *F0t= malloc(sizeof(double)*p*p);
      double *ysq= malloc(sizeof(double)*N*p);
      double *e2= malloc(sizeof(double)*N);
      double *Bsum= malloc(sizeof(double)*p);
      double *sign= malloc(sizeof(double)*p);
      double *Bsdiag= malloc(sizeof(double)*p*p);
      double *gam0t= malloc(sizeof(double)*p*p);
      double *F1t= malloc(sizeof(double)*p*p);
      double *F2= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *s= malloc(sizeof(double)*p);
      double *u= malloc(sizeof(double)*p*p);
      double *ut= malloc(sizeof(double)*p*p);
      double *vtt= malloc(sizeof(double)*p*p);
      double *vttt= malloc(sizeof(double)*p*p);

     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           yy[i*p+j]= y[i+j*N];
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=wt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bs[i*p+j]=invW[i*p+j]*(double)1/sigma0[j];
      for(j=0;j<p;j++){
          for(i=0;i<N;i++){
            wty[j+i*p]=y[j+i*p]*wt[i];
          }
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             xBst[j*N+i]=xBs[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,xBst,&N,&beta,F0t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            ysq[i*p+j]=y[i*p+j]*y[i*p+j];
      for(i=0;i<N;i++){
         e2[i]=0.0;
         for(j=0;j<p;j++){
            e2[i] +=ysq[i*p+j];
         }
         e2[i] = e2[i]*wt[i];
      }
      
      for(j=0;j<p;j++){
         Bsum[j]=0.0;
         for(i=0;i<N;i++){
            Bsum[j] +=Bs[j+i*p]*e2[i];
         }
      }
      for(i=0;i<p;i++){
         for(j=0;j<p;j++){
            Bsdiag[i*p+j]=0.0;
            if(i==j)
              Bsdiag[i*p+i]=Bsum[i];
         }
      }

      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            gam0t[i*p+j]=gam0[i+j*p];
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0t,&p,Bsdiag,&p,&beta,F1t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bsmu[i*p+j]=Bs[i*p+j]*mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]= Bsmu[i*p+j]+alpha0[j]/sigma0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             vt[j*N+i]=v[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,vt,&N,&beta,A0,&p);
      
//      for(i=0;i<p;i++)
//         for(j=0;j<p;j++)
//            F2[i*p+j]=(F0t[i*p+j]-F1t[i*p+j]-A0[i*p+j]);
//      Rprintf("F2 is \n");
//      printmx(F2,p,p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F2[i*p+j]=(F0t[i*p+j]-F1t[i*p+j]-A0[i*p+j])/sumw;
      svd(p,p,F2,s,vtt,u);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            ut[i*p+j]=u[i+j*p];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            vttt[i*p+j]=vtt[i+j*p];
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,ut,&p,vttt,&p,&beta,newgam,&p);
      for(i=0;i<p;i++){
          if(newgam[i*p+i]<0)
            sign[i]=-1.0;
          else
            sign[i]=1.0;
      }
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
           newgam[i*p+j]=newgam[i*p+j]*sign[j];
      if(objgamMS(newgam,y,sigma0,alpha0,mu0,wt,invW,N,p) < objgamMS(gam0,y,sigma0,alpha0,mu0,wt,invW,N,p)){
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            newgam[i*p+j]=gam0[i*p+j];
      }
      free(x); free(vt); free(v); free(Bs); free(wty); free(xBs); free(Bsmu);
      free(xBst); free(F0t); free(ysq); free(e2); free(Bsum); free(sign);
      free(Bsdiag); free(gam0t); free(F1t); free(F2); free(A0); free(s);
      free(u); free(ut); free(vtt);
}
    



     double objgamMS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *w, double *invW, int N, int p){

      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double *F0t= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *x= malloc(sizeof(double)*N*p);
      double *yy= malloc(sizeof(double)*N*p);
      double *v= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *Bs= malloc(sizeof(double)*N*p);
      double *Bsmu= malloc(sizeof(double)*N*p);
      double *xBs= malloc(sizeof(double)*N*p);
      double *xBst= malloc(sizeof(double)*p*N);
      double *wtx= malloc(sizeof(double)*N*p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           yy[i*p+j]= y[i+j*N];

      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      double sumw=0;
      for(i=0;i<N;i++)
         sumw +=w[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            Bs[i*p+j]=invW[i*p+j]*(double)1.0/sigma0[j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
            wtx[j+i*p]=x[j+i*p]*w[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             xBst[j*N+i]=xBs[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,xBst,&N,&beta,F0t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bsmu[i*p+j]=Bs[i*p+j]*mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]=Bsmu[i*p+j]+alpha0[j]/sigma0[j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             vt[j*N+i]=v[j+i*p];
     dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wtx,&p,vt,&N,&beta,A0,&p);
     double val= -1.0*(tr(F0t,p) -2.0*tr(A0,p))/sumw;
     free(F0t); free(A0); free(x); free(v); free(vt); free(Bs);
     free(Bsmu); free(xBs); free(xBst); free(wtx);
     return val;
}
      void updategam1MS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *wt, double *invW, double *newgam, int N, int p){
      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double sumw = 0.0;
      double *x= malloc(sizeof(double)*N*p);
      double *vt= malloc(sizeof(double)*p*N);
      double *v= malloc(sizeof(double)*N*p);
      double *Bsmu= malloc(sizeof(double)*N*p);
      double *Bs= malloc(sizeof(double)*N*p);
      double *xBs= malloc(sizeof(double)*N*p);
      double *xBst= malloc(sizeof(double)*p*N);
      double *wty= malloc(sizeof(double)*N*p);
      double *F0t= malloc(sizeof(double)*p*p);
      double *F1t= malloc(sizeof(double)*p*p);
      double *F1tt= malloc(sizeof(double)*p*p);
      double *F2= malloc(sizeof(double)*p*p);
      double *A0= malloc(sizeof(double)*p*p);
      double *cov1= malloc(sizeof(double)*p*p);
      double *e2= malloc(sizeof(double)*N);
      double *max= malloc(sizeof(double)*N);
      double *row= malloc(sizeof(double)*p);
      double *center= malloc(sizeof(double)*p);
      double *s= malloc(sizeof(double)*p);
      double *sign= malloc(sizeof(double)*p);
      double *u1= malloc(sizeof(double)*N*p);
      double *u= malloc(sizeof(double)*p*p);
      double *ut= malloc(sizeof(double)*p*p);
      double *vtt= malloc(sizeof(double)*p*p);
      double *vttt= malloc(sizeof(double)*p*p);

      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam0,&p,y,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=wt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bs[i*p+j]=invW[i*p+j]*(double)1/sigma0[j];

      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            u1[i*p+j]=y[i*p+j]*invW[i];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
            wty[j+i*p]=y[j+i*p]*wt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            xBs[i*p+j]=x[i*p+j]*Bs[i*p+j];
      for(j=0;j<p;j++)
          for(i=0;i<N;i++)
             xBst[j*N+i]=xBs[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,xBst,&N,&beta,F0t,&p);
      for(i=0;i<N;i++){
         for(j=0;j<p;j++){
             row[j]=Bs[i*p+j];
         }
         max[i]=maxi(row,p);
         e2[i]=max[i]*wt[i];
      }
      for(j=0;j<p;j++)
         center[j]=0.0;
      Covariance(N,p,y,e2,center,cov1);
      
      double e2sum=0.0;
      for(i=0;i<N;i++)
         e2sum +=e2[i];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            cov1[i*p+j]=cov1[i*p+j]*e2sum;
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,gam0,&p,cov1,&p,&beta,F1t,&p);
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             Bsmu[i*p+j]=Bs[i*p+j]*mu0[j];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            v[i*p+j]=Bsmu[i*p+j]+alpha0[j]/sigma0[j];
      for(j=0;j<p;j++)
         for(i=0;i<N;i++)
             vt[j*N+i]=v[j+i*p];
      dgemm_(&notrans,&notrans,&p,&p,&N,&alpha,wty,&p,vt,&N,&beta,A0,&p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F1tt[i*p+j]=F1t[i+j*p];
//      for(i=0;i<p;i++)
//         for(j=0;j<p;j++)
//            F2[i*p+j]=(F0t[i*p+j]-F1tt[i*p+j]-A0[i*p+j]);
//      Rprintf("F2 is \n");
//      printmx(F2,p,p);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            F2[i*p+j]=(F0t[i*p+j]-F1tt[i*p+j]-A0[i*p+j])/sumw;
      svd(p,p,F2,s,vtt,u);
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            ut[i*p+j]=u[i+j*p];
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            vttt[i*p+j]=vtt[i+j*p];
      dgemm_(&notrans,&notrans,&p,&p,&p,&alpha,ut,&p,vttt,&p,&beta,newgam,&p);
      for(i=0;i<p;i++){
          if(newgam[i*p+i]<0.0){
            sign[i]=-1.0;
          }else{
             sign[i]=1.0;}
       }    
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
             newgam[i*p+j]=newgam[i*p+j]*sign[j];
      if(objgamMS(newgam,y,sigma0,alpha0,mu0,wt,invW,N,p) < objgamMS(gam0,y,sigma0,alpha0,mu0,wt,invW,N,p)) {
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            newgam[i*p+j]=gam0[i*p+j];
      }
      free(x); free(vt); free(v); free(Bsmu); free(Bs); free(xBs);
      free(xBst); free(wty); free(F0t); free(F1t); free(F1tt); free(F2);
      free(A0); free(cov1); free(e2); free(max); free(row); free(center);
      free(s); free(sign); free(u); free(ut); free(vtt);
}

      void gig2pMS(double *sr2, int N, int p, double *cpl, double *phi0, double *alpha0, double *W, double *invW, double *logW){
      int i,j;
      double *omega= malloc(sizeof(double)*p);
      double *a1= malloc(sizeof(double)*p);
      double *v1= malloc(sizeof(double)*p);
      double *B1= malloc(sizeof(double)*N*p);
      for(i=0;i<p;i++)
         omega[i]=exp(log(cpl[0+i*2]));
      for(i=0;i<p;i++)
         a1[i]=omega[i]+pow(alpha0[i],2)/phi0[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            B1[i*p+j]=sr2[i*p+j]+omega[j];
      for(i=0;i<p;i++){
          v1[i]=cpl[1+i*2]-((double)1.0/(double)2.0);
      }
      gigp(a1,B1,N,p,v1,W,invW,logW);
      free(omega); free(a1); free(v1); free(B1);
}

       void EMgrstepMSr(double *x, double *mu, double *phi, double *alpha, double **cpl, double **gam, int N, int p, int G, double *pi, int v, int *label, int m){

      int i,g,j;
      double *w= malloc(sizeof(double)*N*G);
      double *mu0= malloc(sizeof(double)*p);
      double *alpha0= malloc(sizeof(double)*p);
      double *phi0= malloc(sizeof(double)*p);
      double *wt= malloc(sizeof(double)*N);
      double *xx= malloc(sizeof(double)*N*p);
     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           xx[i*p+j]= x[i+j*N];
      weightsMS(x,mu,alpha,cpl,gam,phi,N,G,p,v,pi,w);

      if(label)
      combinewk1(w,N,G,label);
      
      for(g=0;g<G;g++){
         for(j=0;j<p;j++){
             mu0[j]=mu[g*p+j];
             phi0[j]=phi[g*p+j];
             alpha0[j]=alpha[g*p+j];
          }
      for(i=0;i<N;i++){
         wt[i]=w[g+i*G];
      }
      updatemaScplpMSr(x,gam[g],mu0,phi0,alpha0,cpl[g],wt,N,p,v,m);
      for(j=0;j<p;j++){
          mu[g*p+j]=mu0[j];
          phi[g*p+j]=phi0[j];
          alpha[g*p+j]=alpha0[j];
      }
   }
      
      for(g=0;g<G;g++){
         pi[g]=0.0;
         for(i=0;i<N;i++){
             pi[g] +=w[g+i*G];
         }
         pi[g] /=N;
      }
      free(w); free(mu0); free(alpha0); free(phi0); free(wt);
}
      





       void updatemaScplpMSr(double *y, double *gam, double *mu0, double *phi0, double *alpha0 , double *cpl, double *wt, int N, int p, int v,int m){
      int i,j;
      char notrans = 'N';
      double alpha = 1.0f;
      double beta = 0.0f;
      double sumw = 0.0;
      double *x= malloc(sizeof(double)*N*p);
      double *yy= malloc(sizeof(double)*N*p);
      double *sr2= malloc(sizeof(double)*N*p);
      double *xmu= malloc(sizeof(double)*N*p);
      double *W= malloc(sizeof(double)*N*p);
      double *invW= malloc(sizeof(double)*N*p);
      double *logW= malloc(sizeof(double)*N*p);
      double *u= malloc(sizeof(double)*N*p);
      double *t= malloc(sizeof(double)*N*p);
      double *newgam= malloc(sizeof(double)*p*p);
      double *mu0new= malloc(sizeof(double)*p);
      double *alpha0new= malloc(sizeof(double)*p);
      double *sigma0new= malloc(sizeof(double)*p);
      double *A= malloc(sizeof(double)*p);
      double *B= malloc(sizeof(double)*p);
      double *C= malloc(sizeof(double)*p);
      double *T= malloc(sizeof(double)*p);
      double *xt= malloc(sizeof(double)*p);
      double *xu= malloc(sizeof(double)*p);
      double *Ax= malloc(sizeof(double)*p);
      double *ax= malloc(sizeof(double)*p);
      double *omega= malloc(sizeof(double)*p);
      double *test1= malloc(sizeof(double)*p*2);
      double *test2= malloc(sizeof(double)*p*3);
      double *row1= malloc(sizeof(double)*2);
      double *row2= malloc(sizeof(double)*3);


     for(i=0;i<N;i++)
        for(j=0;j<p;j++)
           yy[i*p+j]= y[i+j*N];
      if(!wt){
        for(i=0;i<N;i++)
           wt[i]=1;
      }
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,gam,&p,yy,&p,&beta,x,&p);
     for(i=0; i<N;i++){
         for(j=0;j<p;j++){
            x[i*p+j]=pow(x[i*p+j]-mu0[j],2);
            sr2[i*p+j]=x[i*p+j]*1.0/phi0[j];
         }
      }
      gig2pMS(sr2,N,p,cpl,phi0,alpha0,W,invW,logW);
      if(m%2==0)
         updategam2MS(gam,yy,phi0,alpha0,mu0,wt,invW,newgam,N,p);
      else
         updategam1MS(gam,yy,phi0,alpha0,mu0,wt,invW,newgam,N,p);
      dgemm_(&notrans,&notrans,&p,&N,&p,&alpha,newgam,&p,yy,&p,&beta,x,&p);
      for(i=0;i<N;i++)
         sumw +=wt[i];
      for(j=0;j<p;j++){
         A[j]=0.0;
         for(i=0; i<N;i++){
             A[j] +=W[j+i*p]*wt[i];
         }
         A[j] /=sumw;
      }


      for(j=0;j<p;j++){
          B[j]=0.0;
         for(i=0; i<N;i++){
            B[j] +=invW[j+i*p]*wt[i];
         }
         B[j] /=sumw;
      }
      for(j=0;j<p;j++){
          C[j]=0.0;
         for(i=0; i<N;i++){
            C[j] +=logW[j+i*p]*wt[i];
         }
         C[j] /=sumw;
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            u[i*p+j]= (B[j]-invW[i*p+j])*wt[i];
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
            t[i*p+j]=(invW[i*p+j]*A[j]-1.0)*wt[i];
      for(j=0;j<p;j++){
         T[j]=0.0;
         for(i=0;i<N;i++){
            T[j] +=t[j+i*p];
         }
      }
      for(j=0;j<p;j++){
         xt[j]=0.0;
         for(i=0;i<N;i++){
            xt[j] +=x[j+i*p]*t[j+i*p];
         }
      }
      for(j=0;j<p;j++){
         mu0new[j]=xt[j]/T[j];
      }
      for(j=0;j<p;j++){
         xu[j]=0.0;
         for(i=0;i<N;i++){
            xu[j] +=x[j+i*p]*u[j+i*p];
         }
      }
      for(j=0;j<p;j++){
         alpha0new[j]=v*xu[j]/T[j];
      }
      for(i=0;i<N;i++)
         for(j=0;j<p;j++)
             xmu[i*p+j]=pow(x[i*p+j]-mu0new[j],2)*invW[i*p+j];
      for(j=0;j<p;j++){
         Ax[j]=0.0;
         for(i=0;i<N;i++){
            Ax[j] +=xmu[j+i*p]*wt[i];
          }
          Ax[j] /=sumw;
       }
      for(j=0;j<p;j++){
         ax[j]=0.0;
         for(i=0;i<N;i++){
            ax[j] +=x[j+i*p]*wt[i];
         }
         ax[j] /=sumw;
      }
      for(j=0;j<p;j++){
         sigma0new[j]=Ax[j]-2.0*(ax[j]-mu0new[j])*alpha0new[j]+pow(alpha0new[j],2)*A[j];
       }
      for(j=0;j<p;j++)
         omega[j]=exp(log(cpl[0+j*2]));

/*   for(i=0;i<2;i++)*/
       for(j=0;j<p;j++){
           test1[0+j*2]=omega[j];
           test1[1+j*2]=cpl[1+j*2];
       }
       for(j=0;j<p;j++){
          test2[0+j*3]=A[j];
          test2[1+j*3]=B[j];
          test2[2+j*3]=C[j];
       }
       for(j=0;j<p;j++){
          for(i=0;i<2;i++){
             row1[i]=test1[j*2+i];
          }
          for(i=0;i<3;i++){
             row2[i]=test2[j*3+i];
          }
          updateol(row1,row2,2);
          updateol(row1,row2,2);
          for(i=0;i<2;i++){
             cpl[j*2+i]=row1[i];
          }
       }
       for(j=0;j<p;j++){
           if(cpl[1+j*2] <1.0)
             cpl[1+j*2]=1.0;
       }
      for(j=0;j<p;j++)
          phi0[j]=sigma0new[j];
      for(j=0;j<p;j++){
         mu0[j]=mu0new[j];
         alpha0[j]=alpha0new[j];
      }
      for(i=0;i<p;i++)
         for(j=0;j<p;j++)
            gam[i*p+j]=newgam[i*p+j];
   
      free(x); free(sr2); free(xmu); free(W); free(invW); free(logW);
      free(u); free(t); free(newgam); free(mu0new); free(alpha0new);
      free(sigma0new); free(A); free(B); free(C); free(T); free(xt);
      free(xu); free(Ax); free(ax); free(omega); free(test1); 
      free(test2); free(row1); free(row2);

}


void printmx(double *A, int r, int c) {
    int i, j;
        for(i=0; i < r; i++) {
                for(j=0; j < c; j++)
//                        Rprintf("%12.8f ", A[i + r*j]);
                        Rprintf("%12.8f ", A[i*c+j]);
                Rprintf("\n");
        }
        Rprintf("\n");
}

      void combinewk1( double *w, int N, int G, int *label){
      int i,g;
      for(i=0;i<N;i++){
          if(label[i] !=0 && label[i] <= G){
             for(g=0;g<G;g++){
              w[i*G+g]=0.0;}
              w[i*G+label[i]-1]=1.0;
          }else if(label[i] >G){
        for(g=0;g<G;g++)
           w[i*G+g]=0.0;
      }
    }
}           
       















































