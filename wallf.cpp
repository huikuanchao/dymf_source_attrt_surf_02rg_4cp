#include "globals.h"

void wallf(){

int i, j =Dim -1 ;

#pragma omp parallel for private(i)
   for( i = 0; i<nstot ; i++){
	
      if(nstot_flag[i] > 0 ){
	 if(tp[i] != 2){
              if ( x[i][j] < wall_thick )
	         f[i][j] += wall_lamb[tp[i]] ;
	      if ( x[i][j] > (L[j] - wall_thick) )
	         f[i][j] -= wall_lamb[tp[i]]  ;
	  }//tp
	  else{
	   if ( x[i][j] < (wall_thick+Rp) )
	      f[i][j] += wall_lamb[tp[i]] ; 
	   if ( x[i][j] > (L[j] - wall_thick-Rp) )
	     f[i][j] -= wall_lamb[tp[i]]  ;
        

	  }//tp
      
      }
   }

  


}
