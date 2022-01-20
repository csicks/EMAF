cmake_minimum_required(VERSION 3.5.1)
project(EMAF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

add_definitions(-w)

include_directories(align)
include_directories(data)
include_directories(files)
include_directories(utils)

link_libraries(fftw3)
link_libraries(fftw3_threads)
link_libraries(pthread)

add_executable(EMAF
        align/align.cpp
        align/align.h
        align/alignEvaluation.cpp
        align/alignEvaluation.h
        align/polar.cpp
        align/polar.h
        data/array1D.cpp
        data/array1D.h
        data/fileStack.cpp
        data/fileStack.h
        data/image2D.cpp
        data/image2D.h
        data/matrix.cpp
        data/matrix.h
        data/stack3D.cpp
        data/stack3D.h
        files/mrc.cpp
        files/mrc.h
        files/xmippXMD.cpp
        files/xmippXMD.h
        utils/fft.cpp
        utils/fft.h
        utils/logger.cpp
        utils/logger.h
        utils/noise.cpp
        utils/noise.h
        utils/parallel.cpp
        utils/parallel.h
        utils/particle.h
        main.cpp
        )
