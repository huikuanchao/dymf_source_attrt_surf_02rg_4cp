#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "omp.h"
#include "fftw3.h"
using namespace std ;

#define PI   3.141592653589793238462643383

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define min(A,B) ((A)<(B) ? (A) : (B) )
#define Kdelta(i,j) ((i==j) ? 1 : 0 )
#define Dim 3


#ifndef MAIN
extern
#endif
double  eq_phisol,flux_buffer,**avg_rho,wall_lamb[5],wall_thick,lam_L,***W_tsn,***vir_func,***vir_funcpp,***vir_funcpg,**q_noise,Diff_rot,*verlet_a,*verlet_b,p_m,**x_bac,**gn_bac,**grf_bf_x, ***epslon,**real_trq,**trq,**euler_ang,***euler_Q,**euler_adot,quantern4,**euler_q,***euler_B,***euler_A,***dAdphi,***dAdtheta,***dAdpsi,**x, **fnc, **fsol,**f , *tmp, **grid_W, V, L[Dim], dx[Dim], gvol, Lh[Dim], 
     **rho , **w ,  *r_dudr , *tmp2,*tmp3, 
       Diff[5], delt , chiAs,chiBs,chiCs,chiAC,chiBC,chiAB, kappa_p,kappa, C, rho0 , num_averages , 
        Pscalar,Utt , Ubond , Uchi, Ukappa, aver_Ptens[Dim][Dim],Ptens[Dim][Dim], Pvir , phisol,phiHC,phiHA, phiHB,
       *rhoha, *rhow ,*rhosol,*rhohc, *rhohb, *rhoda, *rhodb , *rhot, *rhop, *smrhop,
       **rhosol_t,**rhohc_t,**rhoha_t, **rhohb_t, **rhoda_t, **rhodb_t , **rhop_t,
       CG_ratio, ***Stress_bond_t,*tmp_PP,Stress_PP[Dim][Dim],Stress_Ng[Dim][Dim],Stress_nb[Dim][Dim] , Stress_bonds[Dim][Dim],
       ***sts_buf_pp,***sts_buf_ng,***sts_buf , Range,Range2,Rg, Rg3, phiP, Rp, Xi, Vp, *gammaP,
       *gradwsol[Dim],*gradwC[Dim],*gradwA[Dim], *gradwB[Dim], *gradwP[Dim], 
       *uG, *grad_uG[Dim], *uP, *grad_uP[Dim], *uPG, *grad_uPG[Dim] , mem_use,
       U_chi_gg, U_chi_pg, U_chi_pp, U_kappa_gg, U_kappa_pg, U_kappa_pp,
       sigma, kgraft, **graft_req, *rhoga, **rhoga_t ;

#ifndef MAIN
extern
#endif
int  flux_sp, *nstot_flag,*sol_flag,*nc_flag,*chain_in_buffer,flux_para,wall_para,Nsp ,optm_L, L_dim,L_fren,L_aver,rst_para,uni_sig,nstot, *tp, nA, nB, max_nsol,eq_nsol,nsol,nC ,max_nC,nD, Nha, Nhb,Nhc,Nda, Ndb,
    Nx[Dim], M, nsteps, step, print_freq , nsD,  nsC, nsA, nsB, 
    **grid_inds, pmeorder, spline_weights, lagrange_weights, 
    grid_per_partic, ntypes , stress_freq , buff_size, buff_ind ,
    sample_wait, sample_freq, A_partics , pre_equil_steps ,
    nthreads,
    nsP, nP, Ng, ng_per_partic ;


#ifndef MAIN
extern
#endif
char **xc, tt[80] ;


#ifndef MAIN
extern
#endif
complex<double> *tmp_Ng,**rho_hat,***vir_funcpg_hat ,***vir_funcpp_hat ,***vir_func_hat,*k_wall,*ktmp3,*ktmp2, *ktmp, I, **avg_sk , *grad_uG_hat[Dim] , *grad_uP_hat[Dim], *grad_uPG_hat[Dim] ;

#ifndef MAIN
extern
#endif
long idum;



#ifndef MAIN
extern
#endif
fftw_complex *fin, *fout ;

#ifndef MAIN
extern
#endif
fftw_plan ft_fwd, ft_bck;


#ifndef MAIN
extern
#endif
struct Node {
   int  data ;
   Node* next;
} *nsol_ll_head,*nC_ll_head; 




double ran2(void ) ;
double gasdev2( void ) ;
int cell_stack( int* );
void cell_unstack( int , int* );
void field_gradient( double* , double* , int ) ;
void field_gradient_cdif( double* , double* , int ) ;
void convolve_fields( double*, double*, double* ) ;
void write_grid_data( const char* , double* ) ;
void write_kspace_data( const char* , complex<double>* ) ;

int unstack_stack( int ) ;
void unstack_local( int, int* ) ;
int stack( int* ) ;
int stack_local( int* ) ;
double integrate( double* );
void unstack( int , int* );
double get_k( int , double* ) ;
double get_k_alias( int , double* ) ;
void get_r( int , double* ) ;

void fftw_fwd( double* , complex<double>* );
void fftw_back( complex<double>* , double* );
int remove_dupes( int* , int );

void pbc_vdr( double*, double* , double* );
double pbc_mdr2( double*, double* , double* );
void die(const char *kill);
void matrix_vect( double**, double*,double*);
void matrix_trs_vect(double**,double*,double*);
void wallf();
//Quick routine to kill the program
