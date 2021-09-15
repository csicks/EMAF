cmake_minimum_required(VERSION 3.10)
project(EMAF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_definitions(-w)

include_directories(data)
include_directories(file)
include_directories(utils)
include_directories(evaluation)

link_libraries(fftw3)
link_libraries(fftw3_threads)

add_executable(EMAF
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
        evaluation/alignmentEvaluation.cpp 
        evaluation/alignmentEvaluation.h
        file/mrcFile.cpp
        file/mrcFile.h 
        file/mrcHeader.cpp
        file/mrcHeader.h 
        file/xmippXMD.cpp 
        file/xmippXMD.h
        utils/alignFFT.cpp 
        utils/alignFFT.h
        utils/fft.cpp
        utils/fft.h 
        utils/noise.cpp 
        utils/noise.h
        utils/particle.h 
        utils/polar.cpp 
        utils/polar.h
        main.cpp 
        )
