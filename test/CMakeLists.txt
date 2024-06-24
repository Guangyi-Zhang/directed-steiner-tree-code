cmake_minimum_required(VERSION 3.14)

project(Test LANGUAGES CXX)

# ---- Dependencies ----

CPMAddPackage("gh:doctest/doctest@2.4.9")

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} DST doctest::doctest fmt::fmt)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)

# ---- Add DSTTests ----

enable_testing()

# Note: doctest and similar testing frameworks can automatically configure CMake tests. For other
# testing frameworks add the tests target instead: add_test(NAME ${PROJECT_NAME} COMMAND
# ${PROJECT_NAME})

include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(${PROJECT_NAME})

