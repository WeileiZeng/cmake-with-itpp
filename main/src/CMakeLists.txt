link_directories(${PROJECT_SOURCE_DIR}/lib)
link_libraries(itpp)

include_directories(${PROJECT_SOURCE_DIR}/include/)
add_executable(MainProject main.cpp foo.cpp)
target_link_libraries(MainProject itpp)