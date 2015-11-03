#include "globals.h"
void  add_sol(int );
void add_cchain(int);
int remove_sol();
int find_homoc();
bool deleteNode(struct Node **, Node *);
void addNode( struct Node *, int );

void flux(){

  int tmp_id,tmp_idc[Nhc];
  double cur_ratio = double(nsol)/double(nD*(Nda+Ndb)+nA*Nha + nB*Nhb +nC *Nhc + nsol +nP * ng_per_partic * ( Ng  ) + nP*Vp  );
  if(flux_para <= 0 ) {
	cout<<"wrong flux setup!"<<endl;
       exit(1);
 }

  
  if((flux_para ==  1 )and (cur_ratio < eq_phisol) ){
	tmp_id = find_homoc();
	add_sol( tmp_id);

  }

   if((flux_para ==  2 ) and (cur_ratio > eq_phisol)){
             
           remove_sol();
    
   }


}

int remove_sol(){

    struct Node *cur;
    int j,tmp_rec_nsol[nsol],tmp_nsol=0,i,tp_solind ,ind1, ind2,  ind_t = max_nC*Nhc+nD*(Nda+Ndb) + nA*Nha +nB*Nhb + nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) , ind_t2;

    ind_t2 = ind_t - max_nC*Nhc;
    double tmpz ;

       for(i=0 ;i < nsol ; i++){
	    
	     if(i == 0 ) cur = nsol_ll_head;
	     else cur = cur->next;
             tp_solind = cur->data;
 	     ind1 = ind_t + tp_solind;

             tmpz = x[ind1][Dim-1];//(x[ind1][Dim-1] +x[ind2][Dim-1])/2.0 ;
	     if(tmpz >= flux_buffer){

	       //    deleteNode(&nC_ll_head,cur);	   
		 //  break;
		  tmp_rec_nsol[tmp_nsol] = tp_solind;
		  tmp_nsol +=1;
		  
             }
	     if(tmp_nsol == Nhc ) break;
	}

	if(tmp_nsol < Nhc ){
	 // cout<<"can not find enough solvent in buffer!"<<endl;
	//  exit(1);
	  return 0;
	}


     tmp_nsol = 0;
     for(i=0 ;i < nsol ; i++){
	    
	     if(i == 0 ) cur = nsol_ll_head;
	     else cur = cur->next;
             tp_solind = cur->data;
 	     ind1 = ind_t + tp_solind;

             tmpz = x[ind1][Dim-1];//(x[ind1][Dim-1] +x[ind2][Dim-1])/2.0 ;
	     if(tmpz >= flux_buffer){
                  
	           deleteNode(&nsol_ll_head,cur);	   
		 //  break;
                   tmp_rec_nsol[tmp_nsol] = tp_solind;		  
		  tmp_nsol +=1;
             }
	     if(tmp_nsol == Nhc ) break;
	}



  // delete sol  and ad nchain moleclue
       int tmp_rm_nc = 1;//int(double(tmp_nsol)/double(Nhc));
       
       for(i= 0; i<tmp_rm_nc; i++){

             add_cchain(nC);
             addNode( nC_ll_head, nC);	 

	   for(j =0 ; j<Nhc ; j++){
	    ind1 = tmp_rec_nsol[j];
            sol_flag[ind1] = 0;
	    ind2 =  ind_t + ind1; 
            nstot_flag[ind2] = 0;
	    
	    ind2 = ind_t2+nC*Nhc+j;
	    nstot_flag[ind2] = 1;
	   }

	   nc_flag[nC]  =1 ;
	   nC += 1;
           nsol -= Nhc;
       }
   return 1;
}

void add_cchain(int cind){
	
    struct Node *cur;
    int j,tmp_rec_nsol[nsol],tmp_nsol=0,i,tp_solind ,ind, ind2,  ind_t = max_nC*Nhc+nD*(Nda+Ndb) + nA*Nha +nB*Nhb + nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) , ind_t2;

    ind_t2 = ind_t - max_nC*Nhc;

    ind= ind_t2 + cind*Nhc;

    
    for(j=0 ; j<Dim ; j++){

	if(j < Dim -1)
	   x[ind][j] = ran2() * L[j] ;
	else 
	    x[ind][j] = ran2() *(L[j] - wall_thick -flux_buffer) + flux_buffer;
        
     }
	tp[ ind ] = 3 ;
       
        xc[ind] = "S" ;

	x_bac[ind][0] = -1;
       
	ind +=1;

    for ( i=1 ; i<Nhc ; i++ ) {
        for ( j=0 ; j<Dim ; j++ ) {
	   x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;
       
           if(j!= Dim -1){
	   	if ( x[ind][j] > L[j] )
	      		x[ind][j] -= L[j] ; 
	   	else if ( x[ind][j] < 0.0 )
	      		x[ind][j] += L[j] ;
           }

           if(j==Dim -1){
	   	if ( x[ind][j] < flux_buffer)
		    x[ind][j] += 2.0 * ( flux_buffer - x[ind][j] ) ;
                if ( x[ind][j] > ( L[j] - wall_thick )  )
		  x[ind][j] -= 2.0 * ( x[ind][j] - ( L[j] - wall_thick ) ) ;
           }//j=dim-1
       }

       tp[ ind ] = 3 ;
       xc[ind] = "S" ;
       x_bac[ind][0] = -1;
       ind++ ;
    }

}

int find_homoc(){
       struct Node *cur;
       int i,tp_cind ,ind1, ind2,  ind_t = nD*(Nda+Ndb) + nA*Nha +nB*Nhb + nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) ;     
        double tmpz ;  

       
       for(i=0 ;i < nC ; i++){
	    
	     if(i == 0 ) cur = nC_ll_head;
	     else cur = cur->next;
             tp_cind = cur->data;
 	     ind1 = ind_t + tp_cind*Nhc;
	     ind2 = ind1 + Nhc-1;
             tmpz = (x[ind1][Dim-1] +x[ind2][Dim-1])/2.0 ;
	     if(tmpz >= flux_buffer){

	           deleteNode(&nC_ll_head,cur);	   
		   break;
		   
		   
             }
             if(i == nC-1) {
		cout<<"can not find C chain in the buffer"<<endl;
	        exit(1);
	     }
	}
       return tp_cind;       

}

void  add_sol(int ch_id ){

	
	int i, j,ind, ind_s,ind_t = nD*(Nda+Ndb) + nA*Nha +nB*Nhb + nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) ;

	ind_s = ind_t + max_nC*Nhc +nsol;

        ind = ind_t + ch_id*Nhc;
        for(i = 0 ;i<Nhc ; i++ ){
		for(j =0 ;j<Dim; j++){
                 x[ind_s][j] = x[ind][j];

	        }
	    
	    nstot_flag[ind] = 0;
	    nstot_flag[ind_s] = 1;
            sol_flag[nsol+i] = 1;

	    x_bac[ind_s][0] = -1;
	    tp[ind_s]  = 4 ;
	    xc[ind_s] = "N" ;
	
	    ind_s +=1;
	    ind +=1;
	}
       
        nc_flag[ch_id] = 0;
	nsol += Nhc;
        nC -= 1;
         

}




