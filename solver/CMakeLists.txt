#
# Integrate a SAT solver into the program
#


remove_definitions(-DUSE_CRYPTOMINISAT)
remove_definitions(-DUSE_MINISAT)

if (USE_CRYPTOMINISAT)   # CryptoMiniSat by Mate Soos (http://gitorious.org/cryptominisat)
  add_definitions(-DUSE_CRYPTOMINISAT)
  set(SAT_SOLVER "with CryptoMiniSat 2.9.1 (commit e819ab3236e, date 26/May/2011)"
    PARENT_SCOPE)
  set(sat_LIBS "cryptominisat"
    PARENT_SCOPE)
  set(sat_INCLUDE "cryptominisat"
    PARENT_SCOPE)
  set(sat_MACRO "-DUSE_CRYPTOMINISAT")
  add_subdirectory(cryptominisat)
endif()


if (USE_MINISAT)   # MiniSat by Niklas Een, Niklas Sorensson (http://www.minisat.se/MiniSat.html)
  add_definitions(-DUSE_MINISAT)
  set(SAT_SOLVER "with MiniSat 2.2.0"
    PARENT_SCOPE)
  set(sat_LIBS "minisat"
    PARENT_SCOPE)
  set(sat_INCLUDE "minisat"
    PARENT_SCOPE)
  set(sat_MACRO "-DUSE_MINISAT")
  add_subdirectory(minisat)
endif()

set(sat_DEFINITIONS "${sat_DEFINITIONS} ${sat_MACRO}"
  PARENT_SCOPE)
set(sat_FLAGS "${sat_FLAGS}"
  PARENT_SCOPE)

