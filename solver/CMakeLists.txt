#
# Integrate a SAT solver into the program
#

# We currently support only CryptoMiniSat 2.9.1
# (provisional release, commit 4bd4b69cb81ef, date 1/Apr/2011) by Mate Soos
# http://gitorious.org/cryptominisat

set(SAT_SOLVER "with CryptoMiniSat 2.9.1"
  PARENT_SCOPE)
set(sat_LIBS "cryptominisat"
  PARENT_SCOPE)
set(sat_INCLUDE "cryptominisat"
  PARENT_SCOPE
)
add_subdirectory(cryptominisat)
set(sat_DEFINITIONS "${sat_DEFINITIONS}"
  PARENT_SCOPE)
set(sat_FLAGS "${sat_FLAGS}"
  PARENT_SCOPE)

