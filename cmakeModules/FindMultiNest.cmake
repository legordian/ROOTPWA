#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      cmake module for finding multinest installation
#//      multinest installation location is defined by environment variable $LIBCONFIG
#//
#//      following variables are defined:
#//      MultiNest_DIR              - multinest installation directory
#//      MultiNest_INCLUDE_DIR      - multinest header directory
#//      MultiNest_LIBRARY_DIR      - multinest library directory
#//      MultiNest_LIBS             - multinest library files
#//
#//      Example usage:
#//          find_package(MultiNest 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(MultiNest_FOUND        FALSE)
set(MultiNest_ERROR_REASON "")
set(MultiNest_DEFINITIONS  "")
set(MultiNest_LIBS)


set(MultiNest_DIR $ENV{MULTINEST})
if(NOT MultiNest_DIR)

	set(MultiNest_FOUND TRUE)

	set(_MultiNest_LIB_NAMES "multinest")
	find_library(MultiNest_LIBS
		NAMES ${_MultiNest_LIB_NAMES})
	if(NOT MultiNest_LIBS)
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Cannot find multinest library '${_MultiNest_LIB_NAMES}'.")
	else()
		get_filename_component(MultiNest_DIR ${MultiNest_LIBS} PATH)
	endif()
	unset(_MultiNest_LIB_NAMES)

	set(_MultiNest_HEADER_FILE_NAME "multinest.h")
	find_file(_MultiNest_HEADER_FILE
		NAMES ${_MultiNest_HEADER_FILE_NAME})
	if(NOT _MultiNest_HEADER_FILE)
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Cannot find multinest header file '${_MultiNest_HEADER_FILE_NAME}'.")
	endif()
	unset(_MultiNest_HEADER_FILE_NAME)
	unset(_MultiNest_HEADER_FILE)

	if(NOT MultiNest_FOUND)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} MultiNest not found in system directories (and environment variable MULTINEST is not set).")
	endif()



else()

	set(MultiNest_FOUND TRUE)

	set(MultiNest_INCLUDE_DIR "${MultiNest_DIR}/include")
	if(NOT EXISTS "${MultiNest_INCLUDE_DIR}")
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Directory '${MultiNest_INCLUDE_DIR}' does not exist.")
	endif()

	set(MultiNest_LIBRARY_DIR "${MultiNest_DIR}/lib")
	if(NOT EXISTS "${MultiNest_LIBRARY_DIR}")
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Directory '${MultiNest_LIBRARY_DIR}' does not exist.")
	endif()

	set(_MultiNest_LIB_NAMES "multinest")
	find_library(MultiNest_LIBS
		NAMES ${_MultiNest_LIB_NAMES}
		PATHS ${MultiNest_LIBRARY_DIR}
		NO_DEFAULT_PATH)
	if(NOT MultiNest_LIBS)
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Cannot find multinest library '${_MultiNest_LIB_NAMES}' in '${MultiNest_LIBRARY_DIR}'.")
	endif()
	unset(_MultiNest_LIB_NAMES)

	set(_MultiNest_HEADER_FILE_NAME "multinest.h")
	find_file(_MultiNest_HEADER_FILE
		NAMES ${_MultiNest_HEADER_FILE_NAME}
		PATHS ${MultiNest_INCLUDE_DIR}
		NO_DEFAULT_PATH)
	if(NOT _MultiNest_HEADER_FILE)
		set(MultiNest_FOUND FALSE)
		set(MultiNest_ERROR_REASON "${MultiNest_ERROR_REASON} Cannot find multinest header file '${_MultiNest_HEADER_FILE_NAME}' in '${MultiNest_INCLUDE_DIR}'.")
	endif()
	unset(_MultiNest_HEADER_FILE_NAME)
	unset(_MultiNest_HEADER_FILE)

endif()


# make variables changeable
mark_as_advanced(
	MultiNest_INCLUDE_DIR
	MultiNest_LIBRARY_DIR
	MultiNest_LIBS
	MultiNest_DEFINITIONS
	)


# report result
if(MultiNest_FOUND)
	message(STATUS "Found multinest in '${MultiNest_DIR}'.")
	message(STATUS "Using multinest include directory '${MultiNest_INCLUDE_DIR}'.")
	message(STATUS "Using multinest library '${MultiNest_LIBS}'.")
else()
	if(MultiNest_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested multinest installation:${MultiNest_ERROR_REASON}")
	else()
		if(NOT MultiNest_FIND_QUIETLY)
			message(STATUS "MultiNest was not found:${MultiNest_ERROR_REASON}")
		endif()
	endif()
endif()
