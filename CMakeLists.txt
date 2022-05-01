cmake_minimum_required(VERSION 3.21)
project(calculating_pi_cpp)

set(CMAKE_EXE_LINKER_FLAGS "-static -lgmp")
set(CMAKE_CXX_STANDARD 23)

find_package(Boost 1.60.0 REQUIRED)

#include_directories(C:/boost/include/boost-1_79)
#set(Boost_INCLUDE_DIR /usr/include/boost)
#set(Boost_LIBRARY_DIR c:/mingw/mingw/lib/)
#include_directories(${Boost_INCLUDE_DIR})
find_package(Boost COMPONENTS )



if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(calculating_pi_cpp main.cpp)
    target_link_libraries(calculating_pi_cpp ${Boost_LIBRARIES} -static)
endif()

target_link_libraries(calculating_pi_cpp -static -lgmp)