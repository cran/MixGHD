#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef int __CLPK_integer;
void BesselK(double zr, int nz, double nu, int expon, int nSeq, double *r, double *r1);
void BesselK1(double *zr, int nz, double nu, int expon, int nSeq, double *r, double *r1);
void cov(double *x, double *cov, int N, int p);
double getall(double loglik[], int n);
void rgpar(double *x, int N, int p, int G, double *w, int q);
double randn (double mu, double sigma);
void weightedmean(double *x, int N, int p, double *wt, double *wmean);
void Covariance1(int N, int p, double *x, double *z, double *mu0, double *sampcov);
void svd(int M, int N, double *A, double *s, double *u, double *vtt);
//void svd1(int M, int N, double *A, double *s, double *u, double *vtt)
void eigen(int N, double *A, double *wr, double *vr);
int ginv(__CLPK_integer n, __CLPK_integer lda, double *A, double *B);
void copymx(double *A, int r, int c, int lda, double *C);
void variance(double *x, int N, int p, double *var);
void mahalanobis( int N, int p, double *x, double *mu0, double *cov, double *delta0);
int determinant(double *A, __CLPK_integer k, __CLPK_integer lda, double *res);
void iparM(double *x, int N, int p, double *wt, int q, int G,double *mu0, double *alpha, double *Lambda, double *sigma, double *err, double *cpl);
void rmpar (int p, double *l, int G);
void weights(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double *wg, double *phi, int N, int G, int p,int v, double *pi, double *w);
void ddmsghyp(double *x, double *mu0, double *alpha0, double *phi0, double *cpl0, double *cpl, double *gam, double *wg, double *zlog0, int N, int G, int p);
void ddghyp(double *y, double *mu0, double *alpha0, double *cpl0, double *gam,  double *val1, int N, int p, double *phi0);
void dmsghyp(double *y, double *mu0, double *alpha0, double *phi0, double *cpl, double *gam, double *val2, int N, int p);
void iweights(double *x, double *mu0, double *alpha0, double *phi0, double *cpl, double *cpl0, double *gam, double *wg, int N, int p,int v, double *vg);
void combinewk (double *z, int N, int G, int *labels);
void EMgrstep(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double *wg, double *phi, int N, int G, int p, int v, double *pi, int *label, int m);
void grad(double *x, int N, double *y, double eps, double d, int r, int v, double *vec);
 void updatemaScplp(double *y, double *mu0, double *gam, double *phi0, double *alpha0, int N, int p, double *w, double *cpl, double *cpl0, double *vg, int m);
void updateolU(double *ol, double *ABC, int n);
 void gig2p(double *sr2, int N, int p, double *cpl, double *phi0, double *alpha0,double *W, double *invW, double *logW);
void gigp(double *a, double *B, int N, int p, double *v,double *W, double *invW, double *logW);
void updategam2(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p, double *newgam);
void updategam1(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p, double *newgam);
double objgam(double *gam0, double *y, double *phi0, double *alpha0, double *mu0, double *w, double *vg, double *invW, double *abc, int N, int p);
double llik(double *x, double *mu, double *alpha, double *phi, double *cpl0, double **cpl, double **gam, double *wg, int N, int p, int G, double *pi, int *label);
void MAP(double *x, double *mu, double *alpha, double **cpl, double *cpl0, double **gam, double **sigma, double *wg, double *phi, int N, int G, int p, int v, double *pi, int *label, int *MAP);
double tr(double *A, int N);
int maxi_loc(double *array, int size);
void gig2(double *x, int N, int p, double *mu0, double *alpha0, double *cpl, double *cpl0, double *abc, double *phi0);
void gig(double *b, double a, double v, double *abc, int N);
void mx_vec_mult(int n, int q, double *a, double *b, double *r);
 void vec_mx_mult(int n, int q, double *a, double *b, double *r);
void BAR(double *x, int N, int p, int *l, int G);
double maxi(double *array, int size);
void ddghypFA(double *x, double *mu0, double *alpha0, double *cpl, double *sigma, double *val1, int N, int p,int log1);
void weightsFA(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p,int v, double *pi, double *w);
void gig2FA(double *x, int N, int p, double *mu0, double *alpha0, double *invS, double *cpl1, double *abc);
void gigFA(double *b, double a, double v, double *abc, int N);
void EMgrstepFA(double *x, int N, int p, int G, int v, int q, int *label, double *mu, double *alpha, double *cpl, double **sigma, double **Lambda, double **err, double *pi);
void updatemaScplM(double *x, int N, int p, int q, double *mu0, double *alpha0, double *mu0new, double *alpha0new, double *w, double *cpl1, double *cpl1new, double *sigma, int v, double *Lambda, double *err);
void MAPFA(double *x, double *mu, double *alpha, double *cpl, double **sigma,  int N, int G, int p, int v, double *pi, int *label, int *MAP);
double llikFA(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int p, int G, double *pi);
void rgparGH(double *x, int N, int p, int G, double *w);
void ipar(double *x, int N, int p, double *wt, int G, double *mu0, double *alpha0, double *sigma, double *cpl);
void gig2GH(double *x, int N, int p, double *mu0, double *alpha0, double *invS, double *cpl1, double *abc);
void gigGH(double *b, double a, double v, double *abc, int N);
void updatemaScpl(double *x, int N, int p, double *mu0, double *alpha0,  double *w, double *cpl1, double *sigma, double v, double *mu0new, double *alpha0new, double *cpl1new);
void updateol(double *ol, double *ABC, int n);
void ddghypGH(double *x, double *mu0, double *alpha0, double *cpl, double *sigma, double *val1, int N, int p,int log1);
void weightsGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p,int v, double *pi, double *w);
void EMgrstepGH (double *x, int N, int p, int G, int v, int *label, double *mu, double *alpha, double *cpl, double **sigma, double *pi );
void MAPGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int G, int p, int v, double *pi, int *label, int *MAP);
void dmsghypMS(double *y, double *mu0, double *alpha0, double *cpl, double *gam, double *sigma0, double *val2, int N, int p);
 void weightsMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *phi, int N, int G, int p,int v, double *pi, double *w);
void updatemaScplpMS(double *y, double *gam, double *mu0, double *sigma0, double *alpha0 , double *cpl, double *wt, int N, int p, double *mu0new, double *alpha0new, double *sigma0new, double *newgam, double *cplnew, int v, int m);
double llikMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *sigma, int N, int p, int G, double *pi);
void MAPMS(double *x, double *mu, double *alpha, double **cpl, double **gam, double *sigma, int N, int G, int p, int v, double *pi, int *label, int *MAP);
void updategam2MS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *wt, double *invW, double *newgam, int N, int p);
double objgamMS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *w, double *invW, int N, int p);
void updategam1MS(double *gam0, double *y, double *sigma0, double *alpha0, double *mu0, double *wt, double *invW, double *newgam, int N, int p);
void gig2pMS(double *sr2, int N, int p, double *cpl, double *phi0, double *alpha0, double *W, double *invW, double *logW);
void EMgrstepMSr(double *x, double *mu, double *phi, double *alpha, double **cpl, double **gam, int N, int p, int G, double *pi, int v, int *label, int m);
void updatemaScplpMSr(double *y, double *gam, double *mu0, double *phi0, double *alpha0 , double *cpl, double *wt, int N, int p, int v, int m);
void dgeev_sort(double *Er, double *Ei, double *vr, int N);
void printmx(double *A, int r, int c);
void updatemaScplM1(double *x, int N, int p, int q, double *mu0, double *alpha0, double *mu0new, double *alpha0new, double *w, double *cpl1, double *cpl1new, int v, double *sigma);
void updatemaScplM2(double *x, int N, int p, int q, double *mu0new, double *alpha0new, double *w, double *cpl1new, double *sigma, int v, double *Lambda, double *err);
void Covariance(int N, int p, double *x, double *z, double *mu0, double *sampcov);
 void besselKnuAsym( double x, double nu, int kmax, double *res);
void combinewk1( double *w, int N, int G, int *label);
void EMgrstepMS(double *y, double **gam, double *mu, double *phi, double *alpha, double **cpl, int N, int p, int G, double *pi, int v, int *label, int m);
double llikGH(double *x, double *mu, double *alpha, double *cpl, double **sigma, int N, int p, int G, double *pi);
void gig20(double *sr2sum, int N, int p, double *mu0, double *alpha0, double *cpl, double *cpl0, double *abc, double *phi0);
//extern void dgesdd_( char* jobz, int* m, int* n, double* a,
//                int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
//                               double* work, int* lwork, int* iwork, int* info );





























