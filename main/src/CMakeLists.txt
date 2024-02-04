add_subdirectory(weilei_lib)

link_directories(${PROJECT_SOURCE_DIR}/lib)
link_libraries(itpp)


include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(${PROJECT_SOURCE_DIR}/include/weilei_lib_h)

find_package(OpenMP)


set(MySourcesList "weilei_lib/bp.cpp"  "weilei_lib/dist.cpp" "weilei_lib/mmio.c"  "weilei_lib/lib.cpp"  "weilei_lib/mm_read.cpp"  "weilei_lib/mm_write.cpp"  "weilei_lib/product_lib.cpp")

add_executable(MainProject main.cpp foo.cpp bar.cpp ${MySourcesList})
target_link_libraries(MainProject itpp OpenMP::OpenMP_CXX)

add_executable(examples examples.cpp foo.cpp bar.cpp)
target_link_libraries(examples itpp OpenMP::OpenMP_CXX)


