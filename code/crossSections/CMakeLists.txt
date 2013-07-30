PROJECT(airwaySegmenter)

cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE ( ITK )
IF ( ITK_FOUND )
  INCLUDE(${ITK_USE_FILE})
ELSE ( ITK_FOUND )
  MESSAGE ( FATAL_ERROR "Cannot build without ITK" )
ENDIF ( ITK_FOUND )

FIND_PACKAGE(VTK)
IF (VTK_FOUND)
	INCLUDE (${VTK_USE_FILE})
ELSE(VTK_FOUND)
	MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
ENDIF (VTK_FOUND)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ENDIF(GenerateCLP_FOUND)

#-----------------------------------------------------------------------------
# If choose to build documentation, then search for Doxygen executables.
OPTION(BUILD_DOXYGEN "Build Doxygen Documentation" OFF)
IF(BUILD_DOXYGEN)
  FIND_PACKAGE(Doxygen)
  IF (DOXYGEN_FOUND STREQUAL "NO")
    MESSAGE(FATAL_ERROR "Doxygen not found. Please get a copy http://www.doxygen.org")
  ELSE (DOXYGEN_FOUND STREQUAL "NO")
    # prepare doxygen configuration file
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    ADD_CUSTOM_TARGET(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # cleanup $build/api-doc on "make clean"
    SET_PROPERTY(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES api-doc)

    # add doxygen as dependency to doc-target
    GET_TARGET_PROPERTY(DOC_TARGET doc TYPE)
    IF(NOT DOC_TARGET)
      ADD_CUSTOM_TARGET(doc)
    ENDIF()
    ADD_DEPENDENCIES(doc doxygen)

  ENDIF (DOXYGEN_FOUND STREQUAL "NO")
ENDIF(BUILD_DOXYGEN)

SET (computeAreaAndContourWithMeanNorm_SOURCE computeAreaAndContourWithMeanNorm.cxx)
GENERATECLP(computeAreaAndContourWithMeanNorm_SOURCE computeAreaAndContourWithMeanNorm.xml)

IF(NOT EXECUTABLE_OUTPUT_PATH)
      SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin) 
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)

#INCLUDE_DIRECTORIES( /usr/matlab-R2011a/extern/include )
#LINK_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} /usr/matlab-R2011a/bin/glnxa64 )
#ADD_LIBRARY( SpectralCluster libSpectralCluster.c )

ADD_EXECUTABLE(computeAreaAndContourWithMeanNorm computeAreaAndContourWithMeanNorm.cxx)

TARGET_LINK_LIBRARIES(computeAreaAndContourWithMeanNorm vtkCommon vtkFiltering vtkIO
vtkGenericFiltering )

# vtkRendering vtkHybrid)

# install HTML API documentation and manual pages
#set(DOC_PATH "share/doc/${CPACK_PACKAGE_NAME}-${VERSION}")
#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api-doc/html DESTINATION ${DOC_PATH} )
#
# install man pages into packages, scope is now project root..
#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api-doc/man/man3 DESTINATION share/man/man3/ )
