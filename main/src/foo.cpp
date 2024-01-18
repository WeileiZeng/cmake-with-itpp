#include <iostream>
#include "foo.h"
#include "bar.h"

void foo()
{
  std::cout << "Hello World!\n";
}

void fun1(){
  fun2();
}
