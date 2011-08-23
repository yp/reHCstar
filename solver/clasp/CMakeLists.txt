message("   Building clasp as SAT solver.")

get_property(DEFINITIONS
  DIRECTORY
  PROPERTY COMPILE_DEFINITIONS)
set(sat_DEFINITIONS "${DEFINITIONS}"
  PARENT_SCOPE)
unset(DEFINITIONS)
set(sat_FLAGS ""
  PARENT_SCOPE)

file(GLOB LIBCLASP_SRC *.cpp)

# clasp does not compile under c++0x like the rest of the project
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++98")

include_directories("./")
add_library(clasp STATIC ${LIBCLASP_SRC})
