cmake_minimum_required(VERSION 3.21)
project(Matrix)

set(CMAKE_CXX_STANDARD 20)

set(HEADERS include/Matrix.h include/Vector.h include/QR.h include/EIG.h)

add_executable(Matrix src/TEST_QR.cpp)
