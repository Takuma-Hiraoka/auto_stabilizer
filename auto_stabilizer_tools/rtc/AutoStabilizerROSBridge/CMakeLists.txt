add_compile_options(-std=c++11)
set(CMAKE_BUILD_TYPE Release)

# openrtm has deprecated coding. This causes too many warnings.
#   warning: dynamic exception specifications are deprecated in C++11
add_definitions(-Wno-deprecated)

find_package(catkin REQUIRED COMPONENTS
  rtmbuild
  roscpp
  auto_stabilizer_msgs)

include_directories(${catkin_INCLUDE_DIRS})
rtmbuild_add_library(AutoStabilizerROSBridge SHARED AutoStabilizerROSBridge.cpp)
target_link_libraries(AutoStabilizerROSBridge ${catkin_LIBRARIES})
set_target_properties(AutoStabilizerROSBridge PROPERTIES PREFIX "")
add_dependencies(AutoStabilizerROSBridge RTMBUILD_${PROJECT_NAME}_genrpc)

rtmbuild_add_executable(AutoStabilizerROSBridgeComp AutoStabilizerROSBridgeComp.cpp)
target_link_libraries(AutoStabilizerROSBridgeComp AutoStabilizerROSBridge)
