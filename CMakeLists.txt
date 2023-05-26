cmake_minimum_required(VERSION 3.16.5)
project(NetworkMatcher)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE "Release")

include_directories(NetMatch)

add_executable(NetMatch
        NetMatch/main.cpp
        NetMatch/main.h
        NetMatch/Graph.cpp
        NetMatch/Graph.h
        NetMatch/Bitsprayer.cpp
        NetMatch/Bitsprayer.h
        NetMatch/filesystem.hpp)