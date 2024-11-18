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

# ---- Options ----

OPTION(NO_FAST2_LB "xxx" OFF)
OPTION(NO_FAST3_LB "xxx" OFF)
OPTION(NO_FAST3_PRUNE_U "xxx" OFF)
OPTION(NO_FAST3_PRUNE_SUBTREE "xxx" OFF)

IF(NO_FAST2_LB)
  ADD_DEFINITIONS(-DNO_FAST2_LB)
ENDIF(NO_FAST2_LB)

IF(NO_FAST3_LB)
  ADD_DEFINITIONS(-DNO_FAST3_LB)
ENDIF(NO_FAST3_LB)

IF(NO_FAST3_PRUNE_U)
  ADD_DEFINITIONS(-DNO_FAST3_PRUNE_U)
ENDIF(NO_FAST3_PRUNE_U)

IF(NO_FAST3_PRUNE_SUBTREE)
  ADD_DEFINITIONS(-DNO_FAST3_PRUNE_SUBTREE)
ENDIF(NO_FAST3_PRUNE_SUBTREE)

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

if(USE_EXAMPLE)
  message(STATUS "compiling ./example")
  add_subdirectory(example)
endif()

