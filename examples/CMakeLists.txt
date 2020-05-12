add_executable(1d 1d.cpp)
add_executable(2d 2d.cpp)
add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PUBLIC OpenMP::OpenMP_CXX)
