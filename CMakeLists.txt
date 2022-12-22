cmake_minimum_required(VERSION 3.4)
project(porous_media)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Ofast")
set(PKG_CONFIG_EXECUTABLE "c:/Program Files/pkg-config/pkg-config.exe")

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

find_package(Eigen3 CONFIG REQUIRED)
find_package(jsoncpp CONFIG REQUIRED)

add_executable(porous_media main.cpp)
target_link_libraries(porous_media PRIVATE jsoncpp_lib jsoncpp_object JsonCpp::JsonCpp)
target_link_libraries(porous_media PRIVATE Eigen3::Eigen PkgConfig::gmp)

