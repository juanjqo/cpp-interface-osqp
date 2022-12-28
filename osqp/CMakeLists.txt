cmake_minimum_required(VERSION 3.5)

project(osqp LANGUAGES CXX)
find_package(osqp REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(UNIX AND NOT APPLE)
    FIND_PACKAGE(Eigen3 REQUIRED)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
    ADD_COMPILE_OPTIONS(-Werror=return-type -Wall -Wextra -Wmissing-declarations -Wredundant-decls -Woverloaded-virtual)
endif()

if(APPLE)
    INCLUDE_DIRECTORIES(
        /usr/local/include/
        /usr/local/include/eigen3
        # Most recent versions of brew install here
        /opt/homebrew/include
        /opt/homebrew/include/eigen3)
    ADD_COMPILE_OPTIONS(-Werror=return-type -Wall -Wextra -Wmissing-declarations -Wredundant-decls -Woverloaded-virtual)
    # The library is installed here when using the regular cmake ., make, sudo make install
    LINK_DIRECTORIES(
        /usr/local/lib/)
endif()

set(OSQP_PATH /usr/local/include/osqp)
include_directories(${OSQP_PATH})
add_library(my_osqp /usr/local/lib/libosqp.dylib)

add_executable(osqp main.cpp)
target_link_directories(osqp PRIVATE osqp::osqp)

install(TARGETS osqp
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
