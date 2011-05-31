message("   Building CryptoMiniSat as SAT solver.")

add_definitions(-UUSE_GAUSS)
add_definitions(-DDISABLE_ZLIB)

get_property(DEFINITIONS
  DIRECTORY
  PROPERTY COMPILE_DEFINITIONS)
set(sat_DEFINITIONS "${DEFINITIONS} ${OpenMP_CXX_FLAGS}"
  PARENT_SCOPE)
unset(DEFINITIONS)
set(sat_FLAGS ""
  PARENT_SCOPE)


include_directories("./")

add_library(cryptominisat STATIC
    Solver.cpp
    SolverDebug.cpp
    SolverMisc.cpp
    Gaussian.cpp
    PackedRow.cpp
    XorFinder.cpp
    MatrixFinder.cpp
    VarReplacer.cpp
    ClauseCleaner.cpp
    RestartTypeChooser.cpp
    FailedLitSearcher.cpp
    PartFinder.cpp
    Subsumer.cpp
    PartHandler.cpp
    XorSubsumer.cpp
    StateSaver.cpp
    ClauseAllocator.cpp
    UselessBinRemover.cpp
    OnlyNonLearntBins.cpp
    CompleteDetachReattacher.cpp
    DimacsParser.cpp
    SCCFinder.cpp
    SolverConf.cpp
    ClauseVivifier.cpp
    DataSync.cpp
    BothCache.cpp
)


