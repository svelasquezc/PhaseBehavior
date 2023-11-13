file(GLOB HEADERS CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/*.hpp")

add_library(PhaseBehavior INTERFACE)

# We need this directory, and users of our library will need it too
target_include_directories(PhaseBehavior INTERFACE ../include)

# All users of this library will need at least C++17
target_compile_features(PhaseBehavior INTERFACE cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )