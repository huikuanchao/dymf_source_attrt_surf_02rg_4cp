#include "globals.h"


void initNode(struct Node *head,int n){
	head->data = n;
		head->next =NULL;
}

void addNode( struct Node *head, int n){


	Node *newNode = new Node;

	newNode->data = n;
	newNode->next = NULL;

	Node *cur = head ; 

	while(cur) {
	   if(cur->next == NULL) {
		 cur->next = newNode;
	         return;
	   }
	   cur = cur->next;
	}

}



bool deleteNode(struct Node **head, Node *ptrDel) {
 Node *cur = *head;

 if(ptrDel == *head) {
 	*head = cur->next;
	delete ptrDel;
	return true;

  }

  while(cur) {
	if(cur->next == ptrDel) {
	    cur->next = ptrDel->next;
	    delete ptrDel;
	    return true;

	}
	cur = cur->next;
  }

  cout<<"can not delete Node "<<ptrDel->data<<endl;
  exit(1);
  return false;


}

void initialize_link_list(){
    int i ;

    nC_ll_head = new Node ;
    nsol_ll_head = new Node ;
    initNode(nC_ll_head,0);

   cout<<"here"<<endl;
    for(i= 1; i<nC ;i++){

          addNode( nC_ll_head, i);         

          
    }

    if(nsol > 0 ){

	 initNode(nsol_ll_head,0);
     for(i= 1; i<nsol ;i++){

          addNode( nsol_ll_head, i);
      }

    }

}
