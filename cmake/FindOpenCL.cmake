# Try to find OpenCL
#
# Once done this will define
#
#  OPENCL_FOUND        - system has OpenCL
#  OPENCL_INCLUDE_DIR  - the OpenCL include directory
#  OPENCL_LIB_DIR      - the OpenCL lib directory
#  OPENCL_LIBRARIES    - link these to use OpenCL
#

if ( DEFINED EXTERNAL_PACKAGES_DIR )
	if( WIN32 )
		list( APPEND _opencl_search_dirs "${EXTERNAL_PACKAGES_DIR}/NVidia/CUDA/v8.0/Windows" )
	else( WIN32 )
		list( APPEND _opencl_search_dirs "${EXTERNAL_PACKAGES_DIR}/NVidia/CUDA/v8.0/Linux" )
	endif( WIN32)
endif ()

find_path( OPENCL_INCLUDE_DIR
	NAMES CL/cl.h CL/cl_platform.h
	HINTS ${_opencl_search_dirs}
	PATH_SUFFIXES include
)

find_library( OPENCL_LIBRARIES
	NAMES opencl
	HINTS ${_opencl_search_dirs}
	PATH_SUFFIXES lib lib64 lib32 lib/x64 lib/Win32
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( OpenCL DEFAULT_MSG OPENCL_INCLUDE_DIR OPENCL_LIBRARIES )

mark_as_advanced( OPENCL_INCLUDE_DIR OPENCL_LIBRARIES )
get_filename_component( OPENCL_LIB_DIR ${OPENCL_LIBRARIES} DIRECTORY )
