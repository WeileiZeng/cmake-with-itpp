#include "foo.h"
//#include "bar.h"
#include<cmath>
#include <itpp/itbase.h>

#include <bitset>
//#include "product_lib.h"
#include "weilei_lib.h"
#include <exception>

const int num_cores = 32;
const int debug = 0;
const int n = 12; //n<=30 to avoid negative int for index bx bz

bool descend_col(itpp::GF2mat G){//allow identical cols
  //check if all cols is in descending order, to avoid duplicated case. For rows, call descend_col(G.transpose())
  int c = G.cols();
  int v1=itpp::bin2dec(G.get_col(0));
  int v2;
  for (int i =1; i<c-1;i++){
    v2=itpp::bin2dec(G.get_col(i));
    if (v1<v2)
      return false;
    v1=v2;
  }
  return true;
}

bool increasing_row(itpp::GF2mat G){//exclude identical rows
  //check if all rows is in increasing order, to avoid duplicated case. For columns, call increasing_row(G.transpose())
  int r = G.rows();
  int v1=itpp::bin2dec(G.get_row(0));
  int v2;
  for (int i =1; i<r-1;i++){
    v2=itpp::bin2dec(G.get_row(i));
    if (v1>=v2)
      return false;
    v1=v2;
  }
  return true;
}

bool row_singleton(itpp::GF2mat Gz){
  //check singleton in G: row weight = 1, indicating bas code with distance 1
  itpp::bvec bvec_zero = itpp::zeros_b(Gz.cols());
  bool singleton_Gz=false;
  for ( int i = 0; i < Gz.rows(); i++){
    if ( itpp::BERC::count_errors(bvec_zero, Gz.get_row(i)) == 1){
      return true;
      //      singleton_Gz = true;
      //      break;
      //      std::cout<<"+";
    }
  }
  return false;
}

bool column_zero(itpp::GF2mat Gz){
   //column weight 0 means the single error cannot be detected
   itpp::bvec bvec_zero = itpp::zeros_b(Gz.rows());
// //   bool singleton_Gz=false;
   for ( int i = 0; i < Gz.cols(); i++){
     if ( itpp::BERC::count_errors(bvec_zero, Gz.get_col(i)) == 0){
       return true;
     }
   }
   return false;
}


itpp::GF2mat dec2GF2mat(int dec, int row, int col){
  //convert decimal to matrix with given size.
   int n = row*col;
   itpp::GF2mat M = itpp::GF2mat(itpp::dec2bin( n,dec),false); //convert dec to bin vec
   itpp::GF2mat MM(row,col); //vec to mat
   for ( int i =0;i<row;i++){
     set_submatrix(MM,M.get_submatrix(0,i*col,0,(i+1)*col-1),i,0);
   }
   //	std::cout<<MM<<std::endl;
   return MM;
}

void check_code(int n, itpp::GF2mat & Gx, itpp::GF2mat &Gz){
  CSSCode code;
  code.n = n;
  code.Gx=Gx;
  code.Gz=Gz;
  code.set_up_CxCz();
  code.dist();
  //        int d = code.d;
  if (code.d>=3){
    //          code.k = code.Cx.rows();
    code.k = n-Gx.rows()-Gz.rows();
    std::cout<<code<<std::endl;
    std::cout<<"Gx"<<Gx<<"\nGz"<<Gz<<std::endl;
  }
  return;
}

void iter_J(int r, int rz, int r0, itpp::GF2mat & Gx, itpp::GF2mat & alpha_Gz, itpp::GF2mat &U){
  //run through the rightest block J with size rz * r-rz-r0
  int j_max = std::pow(2,rz * (r-rz-r0));//degree of freedom
  //  std::cout<<"j_max = "<<j_max<<"..."<<std::unitbuf;
  /*
  std::cout<<"j_max = "<<j_max<<" iteration running for"
	   <<" n="<<n
	   <<" rz="<<rz
	   <<" r0="<<r0
	   <<std::endl;*/
  //#pragma omp parallel for schedule(guided) num_threads(num_cores)
  //#pragma omp parallel for num_threads(16)
  //  #pragma omp parallel for schedule(guided) num_threads(16)
  for ( int bj = 1; bj<j_max; bj++){
    itpp::GF2mat J= dec2GF2mat(bj, rz, r-rz-r0);
    /*    
    std::cout<<n<<std::endl;
    std::cout<<rz<<std::endl;
    std::cout<<r0<<std::endl;
    std::cout<<alpha_Gz<<std::endl;
    std::cout<<J<<std::endl;*/
    //reconstruct beta
    itpp::GF2mat omp_alpha_Gz(alpha_Gz);
    set_submatrix(omp_alpha_Gz, J, 0,r0+rz);
    itpp::GF2mat Gz = omp_alpha_Gz*U;
		
    //now construct the code and check distance
    check_code(n,Gx,Gz);
  }
  //  std::cout<<" done"<<std::endl;
  
  return;
}





int test() {
  
  // A CSS code is make of Gx and Gz. For each Gx, run through all possible Gz
    //Gx=(I,M),U=(M^T,I),Gx*U^T=0, Gz=alpha*U


    //number of checks needed to fix any single error
    //    const int r_min =(int) (std::log2(n+1)+0.9999);//for non-degenerate code
    const int r_min = 2; //for degenerate code, e.g. Shor code
    for (int rx = r_min; rx <n-1; rx++){ //number of rows in Gx
      std::cout<<"computing... n="<<n<<", rx="<<rx<<std::endl;
      
      long bx_max = std::pow(2,(n-rx)*rx); //freedom in M
      //      for ( int bx=3803;bx <3803+1;bx++){
      std::cout<<"bx_max="<<bx_max<<std::endl;
#pragma omp parallel for schedule(guided) num_threads(num_cores)
      for ( int bx=1;bx<bx_max;bx++){
	if (debug) std::cout<<"n="<<n 
			    <<",rx="<<rx<<"/"<<n-1
			    <<", bx="<<bx<<"/"<<bx_max
			    << ", bx= " << std::bitset<16>(bx) 
			    << std::endl;
	//print out Gx
	itpp::GF2mat MM = dec2GF2mat(bx, rx, n-rx);

	//	if ( ! is_row_reduced_echelon_form( MM, debug) ) {
	if ( !descend_col(MM) ) {
 	  //std::cout<<"no re\n";
 	  //std::cout<<MM<<std::endl;
	  continue;//skip this case
	}
	else{
	  if (debug) std::cout<<"MM:"<<MM<<std::endl;
	    //	  std::cout<<".";
	}

 	itpp::GF2mat Gx = itpp::gf2dense_eye(rx).concatenate_horizontal(MM);

 	//check singleton in Gx: row weight = 1, incdicating distance = 1; same for column weight 0
 	if (row_singleton(Gx)) continue;
 	else {
	  if (column_zero(Gx)) continue;
	  //	  else std::cout<<"Gx:"<<Gx<<std::endl;
	}
 	itpp::GF2mat U  = MM.transpose().concatenate_horizontal(itpp::gf2dense_eye(n-rx)); //G*U^T=0


 	if (row_singleton(U)) continue;
 	else {
	  if (column_zero(U)) continue;
	  //	  else std::cout<<"U:"<<U<<std::endl;
	}


	//	std::cout<<"bx="<<bx<<"/"<<bx_max<<std::endl;
 	//added for loop
	int r = n-rx; //remained rank
	for ( int rz = r_min; rz < r; rz ++){
	  //add limitation on k
	  int k = n-rz-rx;
	  if (k>1) continue; //only count single logical qubit code for now

	  //	  std::cout<<"rx="<<rx<<",rz="<<rz <<std::endl;

	  //use new method to get Gz


	  //consider zero case
	  itpp::GF2mat alpha_Gz(rz,n-rx);
	  set_submatrix(alpha_Gz, itpp::gf2dense_eye(rz), 0,0); //identity on left
	  iter_J(r, rz, 0,Gx, alpha_Gz, U);
	  
	  //	  if (rx==3 && rz==3 ) std::cout<<"rx="<<rx<<",rz="<<rz <<std::endl;
	  
	  //cases with one or more zero columns
	  for (int r0=1; r0<= r - rz; r0++){   //choose r0 number of zero columns in beta

	    //	    std::cout<<"debug 1: ";
	    
	    itpp::bvec error = itpp::zeros_b(r0+rz);
	    //initialize error
	    for (int i_r0=0;i_r0<r0;i_r0++)
	      error.set(i_r0,1);
	    //std::cout<<error;
	    //error(0)=1;
	    //error.set(r0+rz-1,1);
	    //	    std::cout<<error<<std::endl;;
	    while( true){//all conbinations of zero columns

	      //set corresponding column to zero
	      //set up identity matrix with size rz
	      //std::cout<<" debug 2";
	      if (r0>1) std::cout<<error<<std::endl;

	      itpp::GF2mat alpha_Gz(rz,n-rx);
	      set_submatrix(alpha_Gz, itpp::gf2dense_eye(rz), 0,0); //identity on left
	      int index = rz;//zero column initialized here
	      for (int i = 0;i<r0+rz;i++){ //shift identity matrix by asserting zero columns
		if (error(i)){
		  //permute colmn i and index
		  alpha_Gz.swap_cols(i,index); //only swap? or shift? should be equivalent. need check
		  index ++;
		}
	      }

	      //	      std::cout<<"bx="<<bx<<",";
	      iter_J(r, rz, r0,Gx, alpha_Gz, U);

	      if (!next_error(error, r0+rz,r0))
		break;
	      
	      /*		
	      //run through the rightest block J with size rz * r-rz-r0
	      int j_max = std::pow(2,rz * (r-rz-r0));//degree of freedom
	      for ( int bj = 1; bj<j_max; bj++){
		itpp::GF2mat J= dec2GF2mat(bj, rz, n-rz-r0);
		//reconstruct beta
		set_submatrix(alpha_Gz, J, 0,r0+rz);
		itpp::GF2mat Gz = alpha_Gz*U;

		//now construc t the code and check distance

		//now construct the code
		CSSCode code;
		code.n = n;
		code.Gx=Gx;
		code.Gz=Gz;
		code.set_up_CxCz();
		code.dist();
		//	  int d = code.d;
		if (code.d>2){
		  //	      code.k = code.Cx.rows();
		  code.k = n-rx-rz;
		  std::cout<<code<<std::endl;
		  std::cout<<"Gx"<<Gx<<"\nGz"<<Gz<<std::endl;
		}
	      }
	      */
	    }

	    
	  
	    /* old code
 	  int bz_max = std::pow(2,rz*(n-rx));
	  if (bz_max > 1024*8*4) std::cout<<"rz="<<rz
					  <<", bz_max=2^"<<rz*(n-rx)<<"="<<bz_max<< std::endl;

	  //#pragma omp parallel for schedule(guided) num_threads(num_cores)
#pragma omp parallel for schedule(dynamic, 8) num_threads(num_cores)
	  for ( int bz=1;bz <bz_max-1;bz++){
	    //	  std::cout <<"bz="<<bz<< ", bz= " << std::bitset<24>(bz) << std::endl;
	    itpp::GF2mat alpha_Gz = dec2GF2mat(bz, rz, n-rx);
	    if (!increasing_row(alpha_Gz))
	      continue;

	    //std::cout<<"skip alpha Gz"<<alpha_Gz<<std::endl;	

	    //	  std::cout<<"U"<<U<<std::endl;
	    itpp::GF2mat Gz = alpha_Gz*U;
	    //	  std::cout<<"Gz"<<Gz<<std::endl;
	    //	  Gz = itpp::GF2mat(beta,false)*U;

	    //check singleton in Gz: row weight = 1
	    if (row_singleton(Gz)) continue;
	    if (column_zero(Gz)) {
	      //std::cout<<"." ;
	      continue;
	    }

	    //now construct the code
	    CSSCode code;
	    code.n = n;
	    code.Gx=Gx;
	    code.Gz=Gz;
	    code.set_up_CxCz();
	    code.dist();
	    //	  int d = code.d;
	    if (code.d>2){
	      //	      code.k = code.Cx.rows();
	      code.k = n-rx-rz;
	      std::cout<<code<<std::endl;
	      std::cout<<"Gx"<<Gx<<"\nGz"<<Gz<<std::endl;
	    }
	  //	  std::cout<<code.Gz<<std::endl;
	  //	  std::cout<<code<<std::endl;
	  //	  std::cout<<code<<std::endl;

	  //	  return 0;
	  }
	    */
	  }
	}//end for loop
      
      }
    }
    return 0;
  //  std::cout <<"i="<<i<< ",c = " << std::bitset<16>(i) << std::endl;

  int N = 20000;
  for ( int i =0 ;i<N;i++){
    std::cout <<"i="<<i<< ",c = " << std::bitset<16>(i) << std::endl;
  }
  int i = 0;
  while (i/10000 < 10){
    i++;
    std::cout <<"i="<<i<< ",c = " << std::bitset<24>(i) << std::endl;
  }
  
  return 0;

}

void next_error_test(){
  int n=7;
  itpp::bvec error = itpp::zeros_b(n);
  error.set(0,1);
  error.set(1,1);
  std::cout<<error<<std::endl;;
  while (next_error(error, n , 2)){
    std::cout<<error<<std::endl;
  }
  return;
}


int main()
{
  //  next_error_test(); return 0;
  test();
  return 0;

  
  itpp::GF2mat G(3,3);
  G.set(1,1,1);
  std::cout<<G<<std::endl;
  int n=5;
  itpp::bvec b(n);
  int i=0, max=std::pow(2,n);
  std::cout<<G<<std::endl;
  //  foo();
  

  int x = 0b00001000;
  std::cout<<x<<std::endl;
  std::cout<<x+1<<std::endl;

  test();
  return 0;
}
