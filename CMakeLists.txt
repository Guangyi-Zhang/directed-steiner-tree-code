cmake_minimum_required(VERSION 3.14)

# ---- Project ----

# Note: update this to your new project's name and version
project(
  RootProject
  VERSION 1.0
  LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Add dependencies via CPM ----
# see https://github.com/TheLartians/CPM.cmake for more info
include(cmake/CPM.cmake)

# ---- Subdirectories ----
add_subdirectory(dst)

if(USE_TEST)
  message(STATUS "compiling ./test")
  add_subdirectory(test)
endif()

if(USE_MAIN)
  message(STATUS "compiling ./main")
  add_subdirectory(main)
endif()

