// show examples
#include <itpp/itbase.h>
#include<cmath>

void basic_objects();
void itfile_io();

int main(){
  basic_objects();
  itfile_io();
  return 0;
}

void basic_objects(){

  itpp::GF2mat G(3,5);
  G.set(1,1,1);
  G.set(1,2,1);
  std::cout<<"GF2mat:\n"<<G<<std::endl;
  int n=5;
  itpp::bvec b=itpp::zeros_b(n);
  int i=0, max=std::pow(2,n);
  b.set(2,1);
  std::cout<<"bvec: "<<b<<std::endl;
  std::cout<<"G*b= "<<G*b<<std::endl;
  while ( i < max){
    i++;
    continue;
  }
  for (int i =0; i++;i<5){
    continue;
  }

  return;
}

void itfile_io(){
  //write
  /*
  itpp::GF2mat G(3,5);
  G.set(1,1,1);
  G.set(1,2,1);
  //  std::cout<<"GF2mat:\n"<<G<<std::endl;
  itpp::it_file file("tmp.it");
  file.write_file_header();
  */


    // Declare the it_file class
  itpp::it_file ff;
    // Open a file with the name "it_file_test.it"
    ff.open("it_file_test.it");
    // Create some data to put into the file
    itpp::vec a = itpp::linspace(1, 20, 20);
    // Put the variable a into the file. The Name("a") tells the file class
    // that the next variable shall be named "a".
    ff << itpp::Name("a") << a;
    // Force the file to be written to disc. This is useful when performing
    // iterations and ensures that the information is not stored in any cache
    // memory. In this simple example it is not necessary to flush the file.
    ff.flush();
    // Close the file
    ff.close();
    // Exit program

    //read
    // Declare the it_file class
    itpp::it_file ffr;
    // Open the file "it_file_test.it" for reading
    ffr.open("it_file_test.it");
    // Read the variable a from the file. Put result in vector a.
    itpp::vec ar;
    ffr >> itpp::Name("a") >> ar;
    // Print the result
    std::cout << "a = " << ar << std::endl;
    return ;
}
