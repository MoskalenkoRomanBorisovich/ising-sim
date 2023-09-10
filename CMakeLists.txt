# Works with 3.14 and tested through 3.27
cmake_minimum_required(VERSION 3.14...3.27)

# Project name and a few useful settings. Other commands can pick up the results
project(
  ising
  VERSION 0.1
  DESCRIPTION "Ising Simulation library for c++ and python"
  LANGUAGES CXX)

add_custom_target(bar)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Testing only available if this is the main app
  # Note this needs to be done in the main CMakeLists
  # since it calls enable_testing, which must be in the
  # main CMakeLists.
  include(CTest)
  # enable_testing()
endif()

# The compiled library code is here
add_subdirectory(source)

# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
   AND BUILD_TESTING)
  add_subdirectory(test)
endif()
