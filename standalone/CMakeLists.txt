cmake_minimum_required(VERSION 3.14...3.22)

project(DSTStandalone LANGUAGES CXX)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

CPMAddPackage("gh:Tencent/rapidjson@1.1.0")
CPMAddPackage("gh:gabime/spdlog@1.14.1")

CPMAddPackage(
  GITHUB_REPOSITORY jarro2783/cxxopts
  VERSION 3.2.0
  OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL YES"
)

CPMAddPackage(NAME DST SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create standalone executable ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)

add_executable(${PROJECT_NAME} ${sources})

set(RAPIDJSON_INCLUDE_DIR "${rapidjson_SOURCE_DIR}/include")
set(SPDLOG_INCLUDE_DIR "${spdlog_SOURCE_DIR}/include")
target_include_directories(${PROJECT_NAME} PUBLIC ${RAPIDJSON_INCLUDE_DIR} ${SPDLOG_INCLUDE_DIR})

# set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11 OUTPUT_NAME "DST")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)

target_link_libraries(${PROJECT_NAME} DST::DST cxxopts fmt::fmt)
