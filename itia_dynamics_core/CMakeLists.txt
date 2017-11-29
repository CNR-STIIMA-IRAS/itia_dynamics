cmake_minimum_required(VERSION 2.8.3)
project(itia_dynamics_core)

SET(CMAKE_CXX_FLAGS "-std=c++11 -funroll-loops")
SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3")

find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  kdl_conversions
  itia_futils
  itia_gutils
  itia_mutils
  itia_rutils
  kdl_parser
  std_msgs
  roscpp
)

find_package(orocos_kdl REQUIRED)
find_package(Eigen3 REQUIRED COMPONENTS)
find_package(Boost REQUIRED COMPONENTS system)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES itia_dynamics_core
  CATKIN_DEPENDS roscpp eigen_conversions itia_futils itia_gutils itia_mutils itia_rutils std_msgs kdl_conversions kdl_parser
  DEPENDS system_lib orocos_kdl
)
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
)

add_library(itia_dynamics_core
  src/${PROJECT_NAME}/itia_kin_solver_kdl.cpp
  src/${PROJECT_NAME}/itia_kin_solver_base.cpp
####  src/${PROJECT_NAME}/itia_kin_solver_moveit.cpp
  src/${PROJECT_NAME}/itia_primitives_chain.cpp
  src/${PROJECT_NAME}/itia_primitives_joint.cpp
  src/${PROJECT_NAME}/itia_primitives_link.cpp
)
add_dependencies(itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(itia_dynamics_core ${catkin_LIBRARIES})

add_executable(itia_kin_solver_test src/itia_kin_solver_test.cpp)
add_dependencies(itia_kin_solver_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(itia_kin_solver_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)

add_executable(itia_dyn_solver_test src/itia_dyn_solver_test.cpp)
add_dependencies(itia_dyn_solver_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(itia_dyn_solver_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)

add_executable(itia_inertia_matrix_test src/itia_inertia_matrix.cpp)
add_dependencies(itia_inertia_matrix_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(itia_inertia_matrix_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)

add_executable(itia_dyn_elastic_regressor_test src/itia_dyn_elastic_regressor_test.cpp)
add_dependencies(itia_dyn_elastic_regressor_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(itia_dyn_elastic_regressor_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)


add_executable(jac_speed_test src/jac_speed_test.cpp)
add_dependencies(jac_speed_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(jac_speed_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)


add_executable(vel_speed_test src/vel_speed_test.cpp)
add_dependencies(vel_speed_test itia_dynamics_core ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(vel_speed_test 
                      itia_dynamics_core                       
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${Eigen3_LIBRARIES}
                      ${orocos_kdl_LIBRARIES}
)


