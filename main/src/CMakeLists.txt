link_directories(${PROJECT_SOURCE_DIR}/lib)
link_libraries(itpp)

include_directories(${PROJECT_SOURCE_DIR}/include/)
add_executable(MainProject main.cpp foo.cpp bar.cpp)
target_link_libraries(MainProject itpp)

add_executable(examples examples.cpp foo.cpp bar.cpp)
target_link_libraries(examples itpp)