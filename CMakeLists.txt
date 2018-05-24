cmake_minimum_required( VERSION 3.9 FATAL_ERROR )

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )


project( DiscreteElasticRods )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( EXTERNAL_PACKAGES_DIR "${CMAKE_SOURCE_DIR}/libs/" )

# find Qt components
find_package( Qt5 COMPONENTS Core Widgets OpenGL REQUIRED )
# run moc, uic, rcc automatically when needed
set( CMAKE_AUTOMOC ON )
set( CMAKE_AUTOUIC ON )
set( CMAKE_AUTORCC ON )
# include current src and binary dir to include path to
# properly find files auto-generated from qt tools
set( CMAKE_INCLUDE_CURRENT_DIR ON )

# TODO: CML and DLIB should be submodules
find_package( CML REQUIRED )
find_package( DLIB REQUIRED )

# check for OpenMP support
find_package( OpenMP )
if( OPENMP_FOUND )
	message( STATUS "Found OpenMP support for multithreading" )
else()
	message( WARNING "OpenMP support NOT found. Multithreading might not work" )
endif()

#file(GLOB_RECURSE PROJECT_HEADERS "include/*.h" "include/*.hpp")
#file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp")
#message("PROJECT_HEADERS = ${PROJECT_HEADERS}")
#message("PROJECT_SOURCES = ${PROJECT_SOURCES}")

# headers
set( PROJECT_HEADERS
	include/MainWindow.h
	include/GLWindow.h
	include/GLRenderer.h
	include/GLUtils.h
	include/Camera.h
	include/cml_config.h
	include/config.h
	include/AABB.h
	include/Utils.h
	include/Scene.h
	include/Mesh.h
	include/RenderObject.h
	include/TransformTool.h
	include/ElasticRod.h
	include/Spiral.h
	include/Hair.h
	include/HairGenerator.h
	include/VoxelGrid.h
	include/ObjLoader.h
	include/CollisionShape.h
	include/BasicParser.h
	include/SceneLoader.h
	include/Exporter.h
	include/AnimationBuffer.h
	include/ElasticRodState.h
	include/HairState.h
)

# sources
set( PROJECT_SOURCES
	src/main.cpp
	src/MainWindow.cpp
	src/GLWindow.cpp
	src/GLRenderer.cpp
	src/GLUtils.cpp
	src/Camera.cpp
	src/Utils.cpp
	src/Scene.cpp
	src/Mesh.cpp
	src/RenderObject.cpp
	src/TransformTool.cpp
	src/AABB.cpp
	src/ElasticRod.cpp
	src/Spiral.cpp
	src/Hair.cpp
	src/HairGenerator.cpp
	src/VoxelGrid.cpp
	src/ObjLoader.cpp
	src/CollisionShape.cpp
	src/SceneLoader.cpp
	src/Exporter.cpp
)

add_executable( ${PROJECT_NAME}
	${PROJECT_HEADERS}
	${PROJECT_SOURCES}
)
target_include_directories( ${PROJECT_NAME}
	PRIVATE
		${CMAKE_SOURCE_DIR}/include
		${CML_INCLUDE_DIR}
		${DLIB_INCLUDE_DIR}
)
target_link_libraries( ${PROJECT_NAME}
	Qt5::Core
	Qt5::Widgets
	Qt5::OpenGL
	${CML_LIBRARIES}
	${DLIB_LIBRARIES}
)
set_target_properties( ${PROJECT_NAME}
	PROPERTIES
	AUTOUIC_SEARCH_PATHS "ui ${CMAKE_AUTOUIC_SEARCH_PATHS}"
)

target_compile_definitions( ${PROJECT_NAME}
	PRIVATE
		MULTI_THREAD=1
		#DBUGG=1
)

if( OPENMP_FOUND )
	target_compile_options( ${PROJECT_NAME}
		PRIVATE
			${OpenMP_CXX_FLAGS}
	)
endif()

add_custom_target( COPY_RUNTIME_FILES
	COMMENT "Copying project runtime files to ${PROJECT_BINARY_DIR}"
	COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${PROJECT_BINARY_DIR}/shaders
	COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/assets ${PROJECT_BINARY_DIR}/assets
)
add_dependencies( ${PROJECT_NAME} COPY_RUNTIME_FILES )

# list other files as separate target
file( GLOB_RECURSE PROJECT_UIFORMS "ui/*.ui" )
file( GLOB_RECURSE PROJECT_SHADERS "shaders/*.glsl" )
file( GLOB_RECURSE PROJECT_ASSETS "assets/*" )

add_library( other_files EXCLUDE_FROM_ALL
	${PROJECT_UIFORMS}
	${PROJECT_SHADERS}
	${PROJECT_ASSETS}
)
target_link_libraries( other_files )
set_target_properties( other_files PROPERTIES
	LINKER_LANGUAGE CXX
)
