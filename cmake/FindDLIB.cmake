# Try to find DLIB
#
# Once done this will define
#
#  DLIB_FOUND        - system has DLIB
#  DLIB_INCLUDE_DIR  - the DLIB include directory
#  DLIB_LIB_DIR      - the DLIB lib directory
#  DLIB_LIBRARIES    - link these to use DLIB
#

if ( DEFINED PACKAGE_ROOT_DIR )
	list( APPEND _dlib_search_dirs "${PACKAGE_ROOT_DIR}/dlib-18.9" )
endif ()

find_path( DLIB_INCLUDE_DIR
	NAMES dlib/optimization.h
	HINTS ${_dlib_search_dirs}
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( DLIB DEFAULT_MSG DLIB_INCLUDE_DIR )

set( DLIB_LIB_DIR "" CACHE STRING "DLIB library dir. Using DLIB as header-only library")
set( DLIB_LIBRARIES "" CACHE STRING "DLIB libraries. Using DLIB as header-only library" )

mark_as_advanced( DLIB_INCLUDE_DIR DLIB_LIB_DIR DLIB_LIBRARIES )
