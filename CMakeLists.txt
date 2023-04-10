cmake_minimum_required(VERSION 3.16.5)
project(WeightedSDAforNetworks)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE "Release")

include_directories(BitSprayer)
include_directories(Epi)

add_executable(Testy
        Epi/main.cpp
        Epi/setu.cpp
        Epi/setu.h
        Epi/stat.cpp
        Epi/filesystem.hpp
        Epi/Graph.cpp
        Epi/stat.h
        Epi/Graph.h
        BitSprayer/Bitsprayer.cpp
        BitSprayer/Bitsprayer.h)