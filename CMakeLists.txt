cmake_minimum_required(VERSION 3.10)
project(EMAF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

add_definitions(-w)


set(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(data)
include_directories(file)
include_directories(utils)
include_directories(evaluation)

link_libraries(fftw3)
link_libraries(fftw3_threads)

add_executable(EMAF
        data/array1D.h
        data/image2D.h
        data/stack3D.h
        file/mrcFile.h
        file/mrcHeader.h
        utils/fft.h
        main.cpp utils/fft.cpp data/array1D.cpp file/mrcHeader.cpp file/mrcFile.cpp data/image2D.cpp utils/alignFFT.h utils/alignFFT.cpp utils/polar.h utils/polar.cpp data/stack3D.cpp utils/noise.h utils/particle.h data/matrix.h data/matrix.cpp evaluation/alignmentEvaluation.h data/fileStack.h data/fileStack.cpp file/xmippXMD.h file/xmippXMD.cpp utils/noise.cpp evaluation/alignmentEvaluation.cpp)
