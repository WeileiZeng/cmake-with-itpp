#include <iostream>
#include "foo.h"
#include "bar.h"

void foo()
{
  std::cout << "Hello Foo!\n";
}

void fun1(){
  std::cout << "fun1:I am calling fun2()\n";
  fun2();
}
