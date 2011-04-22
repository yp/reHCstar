message("   Building MiniSAT as SAT solver.")

add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__STDC_FORMAT_MACROS)
add_definitions(-DAVOID_XOR_CLAUSES)

get_property(DEFINITIONS
  DIRECTORY
  PROPERTY COMPILE_DEFINITIONS)
set(sat_DEFINITIONS "${DEFINITIONS}"
  PARENT_SCOPE)
unset(DEFINITIONS)

add_library(minisat STATIC
  Options.cc
  SimpSolver.cc
  Solver.cc
  System.cc
)


