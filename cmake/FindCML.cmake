# Try to find CML
#
# Once done this will define
#
#  CML_FOUND        - system has CML
#  CML_INCLUDE_DIR  - the CML include directory
#  CML_LIB_DIR      - the CML lib directory
#  CML_LIBRARIES    - link these to use CML
#

if ( DEFINED EXTERNAL_PACKAGES_DIR )
	list( APPEND _cml_search_dirs "${EXTERNAL_PACKAGES_DIR}/cml-1_0_2" )
endif ()

find_path( CML_INCLUDE_DIR
	NAMES cml/cml.h cml/version.h
	HINTS ${_cml_search_dirs}
	PATH_SUFFIXES include
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( CML DEFAULT_MSG CML_INCLUDE_DIR )

set( CML_LIB_DIR "" CACHE STRING "CML library dir. CML is header-only library")
set( CML_LIBRARIES "" CACHE STRING "CML libraries. CML is header-only library" )

mark_as_advanced( CML_INCLUDE_DIR CML_LIB_DIR CML_LIBRARIES )
