cmake_minimum_required(VERSION 3.4)
project(porous_media)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Ofast -march=native")

find_package(Eigen3 CONFIG REQUIRED)
find_package(jsoncpp CONFIG REQUIRED)

add_executable(porous_media main.cpp)
target_link_libraries(porous_media PRIVATE jsoncpp_lib jsoncpp_object JsonCpp::JsonCpp)
target_link_libraries(porous_media PRIVATE Eigen3::Eigen)

add_executable(dynamic_porous_media dynamic_main.cpp)
target_link_libraries(dynamic_porous_media PRIVATE jsoncpp_lib jsoncpp_object JsonCpp::JsonCpp)