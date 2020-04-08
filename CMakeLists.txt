cmake_minimum_required (VERSION 3.9)

add_compile_options(-std=c++11)
add_compile_options(-O3)

set(ROOT ${PROJECT_SOURCE_DIR})
set(SRC "${ROOT}/src")
include_directories("${SRC}")

find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(QSpline ${SRC}/quinticSpline.cpp)

add_executable(example00 ${ROOT}/examples/00.cpp)
target_link_libraries(example00 PUBLIC QSpline
${EIGEN3_LIBS}
)

add_executable(example01 ${ROOT}/examples/01.cpp)
target_link_libraries(example01 PUBLIC QSpline
${EIGEN3_LIBS}
)