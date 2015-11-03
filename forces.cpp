#include "globals.h"
void charge_grid( void ) ;
void bonds( void ) ;


void forces() {

  int i,j, m, gind, t1, t2 ;

  charge_grid() ;


  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////

#pragma omp parallel for private(j)
  for ( i=0 ; i<nstot ; i++ )
    for ( j=0 ; j<Dim ; j++ )
      f[i][j] = 0.0 ;
 
#pragma omp parallel for private(j)
  for ( i=0 ; i<M ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwsol[j][i]=gradwC[j][i] = gradwA[j][i] = gradwB[j][i] = gradwP[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////
  fftw_fwd( rho[0] , ktmp ) ;
  if ( chiAB != 0.0 ) {
    // A acting on B //
 
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwB[j][i] += tmp[i] * chiAB / rho0 ;
    }
 
   }//chiAB

 // A acting on C // 
   if ( chiAC != 0.0 ) {
 
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwC[j][i] += tmp[i] * chiAC / rho0 ;
 

    }
   }

   // A acting on sol // 
   if(chiAs != 0.0 and nsol >0 ){
     for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwsol[j][i] += tmp[i] * chiAs / rho0 ;
    }
 
   }//chiAs
  
 // A Monomers acting on particles //
  if(nP > 0){
   for ( j=0 ; j<Dim ; j++ ) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;

    }
  }//nP

  //end of A acting part //

    // B acting on A //
    fftw_fwd( rho[1] , ktmp ) ;
   
   if ( chiAB != 0.0 ) {
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwA[j][i] += tmp[i] * chiAB / rho0 ;
    }
  }//chiAB
 
 //B acting on C//
  if(chiBC != 0.0 ){
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwC[j][i] += tmp[i] * chiBC / rho0 ;
    }
 
  }//chiBC

//B acting on sol
  if(chiBs != 0.0 and nsol > 0){
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwsol[j][i] += tmp[i] * chiBs / rho0 ;
    }

  }//chiBs

// B Monomers acting on particles //    
  if(nP >0){
    for ( j=0 ; j<Dim ; j++) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;

    }
  }//nP

//end of B acting part//

fftw_fwd( rho[3] , ktmp ) ;

   if ( chiAC != 0.0 ) {
 
 
    // C acting on A //
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwA[j][i] += tmp[i] * chiAC / rho0 ;
    }
  }
  
   if ( chiBC != 0.0 ) {
 
 
    // C acting on B //
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwB[j][i] += tmp[i] * chiBC / rho0 ;
    }
  }
  
  //C acting on sol//

  if (chiCs != 0.0 and nsol >0 ){
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwsol[j][i] += tmp[i] * chiCs / rho0 ;
    }
 

  }//chiCs
  
  // C Monomers acting on particles //
  if(nP >0){
    for ( j=0 ; j<Dim ; j++) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAC : 0.0 ) ) / rho0 ;

    }
  }//nP
 

//end of C acting part//


  fftw_fwd( rho[4] , ktmp ) ;

   if ( chiAs != 0.0 and nsol > 0 ) {
 
 
    // sol acting on A //
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwA[j][i] += tmp[i] * chiAs / rho0 ;
    }
  }

   if ( chiBs != 0.0 and nsol > 0 ) {
 
 
    // sol acting on B //
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwB[j][i] += tmp[i] * chiBs / rho0 ;
    }
  }

   if ( chiCs != 0 and nsol > 0 ){
    // sol acting on C //
    
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;
 
      fftw_back( ktmp2 , tmp ) ;
 
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwC[j][i] += tmp[i] * chiCs / rho0 ;
    }
  }

 
 
 // sol  Monomers acting on particles //
 if(nP > 0){
  for ( j=0 ; j<Dim ; j++) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAs : 0.0 ) ) / rho0 ;

    }
  }//nP

  //end of sol acting part//

  // Compressibility contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] = rhot[i] ;

  fftw_fwd( tmp , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ ) {
      gradwA[j][i] += tmp[i] * kappa / rho0 ;
      gradwB[j][i] += tmp[i] * kappa / rho0 ;
      gradwC[j][i] += tmp[i] * kappa / rho0 ;
      gradwsol[j][i] += tmp[i] * kappa / rho0 ;
    }

  }


  ///////////////////////////////////////////
  // Accumulate the particle contributions //
  ///////////////////////////////////////////
  if ( nP > 0 ) {

    // Particle-particle //
    fftw_fwd( rho[2] , ktmp ) ;


    if(nP >1){
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uP_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * kappa_p / rho0 ;
    }
    }//nP>1
    // Particles acting on monomers //
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ ) {
        gradwA[j][i] += tmp[i] * kappa / rho0 ;
        gradwB[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;
        gradwC[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAC : 0.0 ) ) / rho0 ;
       gradwsol[j][i] += tmp[i] * (kappa + ( A_partics ? chiAs : 0.0 ) )  / rho0 ;
      }
    }


  }// if ( nP > 0 )




  //////////////////////////////////////////////////////
  // Accumulate the nonbonded forces on each particle //
  //////////////////////////////////////////////////////
#pragma omp parallel for private(m,j)
  for ( i=0 ; i<nstot ; i++ ) {

    if( nstot_flag[i] >0 ){
    for ( m=0 ; m < grid_per_partic ; m++ ) {
      
      gind = grid_inds[ i ][ m ] ;

      for ( j=0 ; j<Dim ; j++ ) {
        if ( tp[i] == 0 )
          f[i][j] -= gradwA[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 1 )
          f[i][j] -= gradwB[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 2 )
          f[i][j] -= gradwP[ j ][ gind ] * grid_W[i][m] ;
      	else if (  tp[i] ==3 )
	  f[i][j] -= gradwC[ j ][ gind ] * grid_W[i][m] ;
        else if (  tp[i] ==4 )
	  f[i][j] -= gradwsol[ j ][ gind ] * grid_W[i][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[i][j] *= gvol ;
   }//if nstot
   else{

   }
  }


  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  bonds() ;
  if(wall_para > 0)  wallf();
  
}
