#include "foo.h"
//#include "bar.h"
#include<cmath>
#include <itpp/itbase.h>

int main()
{
  itpp::GF2mat G(3,3);
  G.set(1,1,1);
  std::cout<<G<<std::endl;
  int n=5;
  itpp::bvec b(n);
  int i=0, max=std::pow(2,n);
  std::cout<<G<<std::endl;
  while ( i < max){
    i++;
    continue;
  }
  for (int i =0; i++;i<5){
    continue;
  }
  foo();
  return 0;
}
