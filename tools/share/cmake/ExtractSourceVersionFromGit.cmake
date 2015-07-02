# - ExtractSourceVersionFromGit
# This module extract the commit-id of the last commit of the current
# git repository. This code sets the following variables:
#  
#  APPLICATION_SOURCE_VERSION:  commit-id of the last commit of the repository
#

#
# Get the source information
#
set (APPLICATION_SOURCE_VERSION "unknown")

execute_process(COMMAND sh -c "sh thirdparty/autorevision.sh -t sh -o VERSION; echo echo \\\$VCS_TAG-\\\$VCS_TICK-g\\\$VCS_SHORT_HASH"
  COMMAND sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE APPLICATION_SOURCE_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if (${APPLICATION_SOURCE_VERSION} STREQUAL "unknown")
message (STATUS "!!! The source code is NOT in a git repository.")
else (${APPLICATION_SOURCE_VERSION} STREQUAL "unknown")
message (STATUS "Application source code version: " ${APPLICATION_SOURCE_VERSION})
endif (${APPLICATION_SOURCE_VERSION} STREQUAL "unknown")
