#add_subdirectory(weilei)

link_directories(${PROJECT_SOURCE_DIR}/lib)
link_libraries(itpp)

include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(${PROJECT_SOURCE_DIR}/include/weilei_lib_h)
#add_executable(MainProject main.cpp foo.cpp bar.cpp)
#target_link_libraries(MainProject itpp)

#add_executable(examples examples.cpp foo.cpp bar.cpp)
#target_link_libraries(examples itpp)


add_executable(test_lib test_lib.cpp   bp.cpp  dist.cpp mmio.c  lib.cpp  mm_read.cpp  mm_write.cpp  product_lib.cpp)
target_link_libraries(test_lib itpp)