#include "foo.h"
//#include "bar.h"
#include<cmath>
#include <itpp/itbase.h>

#include <bitset>
//#include "product_lib.h"
#include "weilei_lib.h"
#include <exception>

const int num_cores = 5;
const int debug = 0;
const int n = 20; //n<=30 to avoid negative int for index bx bz

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
      //      std::cout<<bx_max<<std::endl;
      for ( int bx=1;bx <bx_max;bx++){
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

 	//added for loop
	for ( int rz = r_min; rz < n-rx; rz ++){
	  //add limitation on k
	  int k = n-rz-rz;
	  if (k>1) continue;


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

int main()
{
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
