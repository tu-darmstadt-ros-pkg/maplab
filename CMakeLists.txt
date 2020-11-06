cmake_minimum_required(VERSION 2.8.3)
project(imu_integrator_rk4)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

cs_add_library(${PROJECT_NAME}
  src/imu-integrator.cc
  src/imu-integrator-eigen.cc
)
target_link_libraries(${PROJECT_NAME} pthread)

cs_install()
cs_export()
