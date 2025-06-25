CMAKE_MINIMUM_REQUIRED(VERSION 3.5...3.26)
PROJECT(dqrobotics-interface-osqp)
set (CMAKE_CXX_STANDARD 11)

################################################################
# INSTALL HEADERS
################################################################

INSTALL(FILES
    include/dqrobotics/solvers/DQ_OSQPSolver.h
    DESTINATION "include/dqrobotics/solvers")
