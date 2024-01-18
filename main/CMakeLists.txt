# cmake_minimum_required(VERSION 2.4)
cmake_minimum_required(VERSION 2.8.6)
# cmake_minimum_required(VERSION 3.8)

project(hello_itpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_subdirectory(src)



#include_directories(${PROJECT_SOURCE_DIR}/include)
#link_directories(${PROJECT_SOURCE_DIR}/lib)
#link_libraries(itpp)
#add_executable(app src/main.cpp src/foo.cpp)  # be sure there's exactly one main() function in the source files
#target_link_libraries(app itpp )



enable_testing()
add_test(NAME HelloWorld
         COMMAND ./bin/MainProject)

add_custom_target(TestRun ./bin/MainProject)
add_custom_command(OUTPUT null.txt COMMAND ./bin/MainProject)