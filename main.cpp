#define MAIN
#include "globals.h"
void flux();
void write_film_gro();
void adj_L(void);
void update_positions( void ) ;
void update_euler(void);
void initialize( void ) ;
void write_np(void) ;
void write_gro( void ) ;
void write_rst_gro(void );
void write_quaternions(void) ;
void write_grid( void ) ;
void forces( void ) ;
void torque(void);
double integrate( double* ) ;
void write_stress( void ) ;
void bond_stress( void ) ;
void calc_Unb( void ) ;
void calc_stress();
void read_input() ;

int main( int argc , char** argv ) {

  int L_flag=0,i,j,k,l;

  if ( argc > 1 and !strcmp( "-nt" , argv[1] ) ){
    nthreads = atoi( argv[2] ) ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }
  else {
    nthreads = 1 ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
   }

  read_input() ;

   if(argc == 5 ) {
     flux_sp =  atoi(argv[4]);
     cout<<"new flus sp "<<flux_sp<<endl;
   }
  
  initialize() ;

  write_gro( ) ;
  write_quaternions( ) ;
   write_film_gro();
  // Save chi for pre-equilibration steps //
  double tmp_ang,chi_bkp = chiAB ;

  FILE *otp ,*otpL;
  otp = fopen( "data.dat" , "w" ) ;
  otpL = fopen("box_L.dat","w");

  printf("Entering main loop!\n") ; fflush( stdout ) ;
  for ( step = 0 ; step < nsteps ; step++ )  {

    if ( step < pre_equil_steps )
      chiAB = 0.0 ;
    else
      chiAB = chi_bkp ;
   

    forces() ;
  

    if(sigma>0){
       torque();
    }

   // if(step >0 || rst_para == 1)
   	update_positions() ;
   // else
    //	update_positions_init() ;
    
    if(sigma>0){
    update_euler();
    }
   
    if( ( flux_para == 1) and (step % flux_sp == 0)){


	 flux();

    }
    else if((flux_para == 2 ) and (step % flux_sp == 0) and (max_nC> nC ) and (nsol>Nhc)){

         flux();	
    }
    else{

    }


	 

    if ( stress_freq > 0 && step % stress_freq == 0 ) {
      calc_stress() ;

      for ( j=0 ; j<Dim ; j++ ){ 
        for ( k=0 ; k<Dim ; k++ ) {
          sts_buf[buff_ind][j][k]= Rg3*Ptens[j][k];//( j<Dim ? Stress_bonds[j][k]:0.0) ;
           sts_buf_pp[buff_ind][j][k] = Rg3*Stress_PP[j][k];
      	   sts_buf_ng[buff_ind][j][k] = Rg3*Stress_Ng[j][k];
            
	    if( ((L_fren-L_aver) < L_flag)  and (optm_L >0) and (j ==k))
	        aver_Ptens[j][j] += Rg3*Ptens[j][j];


	/*for ( k=0 ; k<nP;  k++ ){
		//euler_adot[k][j] = euler_q[k][j];	
		sts_buf[buff_ind][j][k+Dim] = euler_q[k][j];
	}*/
	}
      }
      if(((L_fren-L_aver) < L_flag) and   (optm_L >0))
      aver_Ptens[0][1] +=1;
 
      buff_ind++ ;
    }

   if(optm_L>0  ){
    	if( step > pre_equil_steps )
       		L_flag +=1 ;

    	if(L_flag == L_fren){
       		adj_L();
		L_flag = 0 ;
               
		
		
	        fprintf( otpL , "%d %lf %lf %lf \n" ,step ,L[0],L[1],L[2]);fflush( otpL ) ;
	}		
   }//optm_L

    if(step % sample_freq == 0){
	write_np();
        //write_film_gro();
        char nm[20];

	if(nsol > 0){
	   sprintf( nm , "./frame/rhosol.frame%d.dat" , step ) ;
	   write_grid_data( nm , rhosol ) ;
	}
	if(nD > 0){
	   sprintf( nm , "./frame/rhoda.frame%d.dat" , step ) ;
	    write_grid_data( nm , rhoda ) ;
	   sprintf( nm , "./frame/rhodb.frame%d.dat" , step ) ;
	               write_grid_data( nm , rhodb ) ;
	}
    }  
    
    if ( step > sample_wait && step % sample_freq == 0 ) {
     /* fftw_fwd( rho[0] , ktmp ) ;
      for ( i=0 ; i<M ; i++ ) {
        avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
      }

      if ( nP > 0 ) {

        fftw_fwd( rho[2] , ktmp ) ;
        for ( i=0 ; i<M ; i++ ) {
          avg_sk[2][i] += ktmp[i] * conj( ktmp[i] ) ;
        }
      }*/
     /* for ( i=0 ; i<M ; i++ ) {
        avg_rho[0][i] += rho[0][i];
	avg_rho[1][i] += rho[1][i];
        //avg_rho[3][i] += rho[3][i];
       
      }*/
	


      num_averages += 1.0 ;
    }


    if ( step % print_freq == 0 || step == nsteps-1 ) {
      printf("step %d of %d  Ubond: %lf\n" , step , nsteps , Ubond ) ;
      fflush( stdout ) ;
      write_gro() ;
      write_rst_gro();
      write_quaternions();


      if ( stress_freq > 0 )
        write_stress() ;


      write_grid_data( "rhoda.dat" , rhoda ) ;
      write_grid_data( "rhodb.dat" , rhodb ) ;
      
      if ( nA > 0.0 )
        write_grid_data( "rhoha.dat" , rhoha ) ;

      if ( nB > 0.0 )
        write_grid_data( "rhohb.dat" , rhohb ) ;
      if(  nC > 0.0 )
         write_grid_data( "rhohc.dat" , rhohb ) ;

      if ( nP > 0.0 ) 
        write_grid_data( "rhop.dat" , rhop ) ;

      if ( step > sample_wait ) {
       /* for ( i=0 ; i<M ; i++ ) 
          ktmp2[i] = avg_sk[0][i] / num_averages ;
        write_kspace_data( "avg_sk_A.dat" , ktmp2 ) ;
        
        if ( nP > 0 ) {
          for ( i=0 ; i<M ; i++ )
            ktmp2[i] = avg_sk[2][i] / num_averages ;
          write_kspace_data( "avg_sk_np.dat" , ktmp2 ) ;
        }*/
        /*	
        for ( i=0 ; i<M ; i++ )
	     tmp[i] = avg_rho[0][i]/ num_averages ;
         write_grid_data("avg_typeA.dat",tmp);	

	for ( i=0 ; i<M ; i++ )
	     tmp[i] = avg_rho[1][i]/ num_averages ;
         write_grid_data("avg_typeB.dat",tmp);	

        for ( i=0 ; i<M ; i++ )
	     tmp[i] = avg_rho[3][i]/ num_averages ;
         write_grid_data("avg_typeC.dat",tmp);	
	*/

      }

      calc_Unb() ;

      fprintf( otp , "%d %lf %lf %lf %lf %lf %lf %lf\n" , step , Ubond , U_chi_gg, U_kappa_gg,U_chi_pg,U_kappa_pg,U_kappa_pp , Utt) ;
      fflush( otp ) ;

    }// if step % print_Freq == 0

  }

  fclose( otp ) ;

  return 0 ;

}
