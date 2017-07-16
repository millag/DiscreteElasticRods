cmake_minimum_required( VERSION 3.1 )

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )


project( DiscreteElasticRods )

set( PACKAGE_ROOT_DIR "${CMAKE_SOURCE_DIR}/libs/" )

# run moc automatically when needed.
set( CMAKE_AUTOMOC ON )
# find the QtWidgets library
find_package( Qt5Core )
find_package( Qt5Widgets )
find_package( Qt5OpenGL )

# TODO: CML and DLIB should be submodules
find_package( CML REQUIRED )
find_package( DLIB REQUIRED )

set( PROJECT_HEADERS
	include/MainWindow.h
	include/GLWindow.h
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

set( PROJECT_SOURCES
	src/main.cpp
	src/MainWindow.cpp
	src/GLWindow.cpp
	src/GLUtils.cpp
	src/Camera.cpp
	src/Utils.cpp
	src/Scene.cpp
	src/Mesh.cpp
	src/RenderObject.cpp
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

#file(GLOB_RECURSE PROJECT_HEADERS "include/*.h" "include/*.hpp")
#file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp")
#message("PROJECT_HEADERS = ${PROJECT_HEADERS}")
#message("PROJECT_SOURCES = ${PROJECT_SOURCES}")

add_executable( ${PROJECT_NAME}
	${PROJECT_HEADERS} ${PROJECT_SOURCES}
	)
target_include_directories( ${PROJECT_NAME}
	PRIVATE
		. ${CMAKE_SOURCE_DIR}/include ${CML_INCLUDE_DIR} ${DLIB_INCLUDE_DIR}
	)
target_link_libraries( ${PROJECT_NAME}
	Qt5::Core Qt5::Widgets Qt5::OpenGL ${CML_LIBRARIES} ${DLIB_LIBRARIES} )

add_custom_target( COPY_RUNTIME_FILES
	COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${PROJECT_BINARY_DIR}/shaders
	COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/assets ${PROJECT_BINARY_DIR}/assets
	COMMENT "Copying project runtime files to ${PROJECT_BINARY_DIR}"
	)
add_dependencies( ${PROJECT_NAME} COPY_RUNTIME_FILES )
