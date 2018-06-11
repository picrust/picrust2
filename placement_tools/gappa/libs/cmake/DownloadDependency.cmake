# gappa - Genesis Applications for Phylogenetic Placement Analysis
# Copyright (C) 2017-2018 Lucas Czech and HITS gGmbH
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@h-its.org>
# Exelixis Lab, Heidelberg Institute for Theoretical Studies
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# --------------------------------------------------------------------------------------------------
#   Purpose
# --------------------------------------------------------------------------------------------------

# This is a helper script that downloads dependencies if they are not there.
# This happens if the REPOSITORY was cloned without --recursive, or directly downloaded from github.

# --------------------------------------------------------------------------------------------------
#   Settings
# --------------------------------------------------------------------------------------------------

# Path to the helper script for doing the actual download. Has to be set outside of the functions,
# see https://stackoverflow.com/a/12854575/4184258
SET(DOWNLOAD_EXTERNAL_FILE_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/DownloadExternalFile.cmake)

# --------------------------------------------------------------------------------------------------
#   Download Functions
# --------------------------------------------------------------------------------------------------

# Download a dependency if it is not found.
#
# The function expects four parameters:
#
#  - LIBPATH  : Path to the libracy dir where dependencies are stored.
#  - LIBNAME  : Name of the dependency, that is, the name of its main directory within the ${LIBPATH}.
#  - TESTFILE : A testfile to check if the dependency is already there.
#  - LIBURL   : The URL to download from if it is not there.
#
# For checking if it is already there, we use a test-file that we know should exist if the
# dependency is already there. The file is search for in ${LIBPATH}/${LIBNAME}/${TESTFILE}.
function( DOWNLOAD_DEPENDENCY LIBPATH LIBNAME TESTFILE LIBURL )

    IF( NOT EXISTS ${LIBPATH}/${LIBNAME}/${TESTFILE} )
        message (STATUS "${LIBNAME} not found")
        message (STATUS "${ColorBlue}Downloading ${LIBNAME} from ${LIBURL}${ColorEnd}")

        # If the file was not found, we download and unpack it (at configure time). This roughly follows
        # https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

        # The DownloadExternalFile.cmake script contains two variables to be replaced, which we set here.
        # In order to not replace any other variables used by the script itself,
        # we only replace @-variables, see https://cmake.org/cmake/help/v3.0/command/configure_file.html
        # Then, we use the configure_file command to make a copy of the script
        # that we then execute to download the dependency.
        SET(DEPENDENCY_PATH ${LIBPATH})
        SET(DEPENDENCY_NAME ${LIBNAME})
        SET(DEPENDENCY_URL  ${LIBURL})
        configure_file(
            ${DOWNLOAD_EXTERNAL_FILE_SCRIPT}
            ${CMAKE_BINARY_DIR}/${LIBNAME}Download/CMakeLists.txt
            @ONLY
        )

        execute_process( COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${LIBNAME}Download
        )

        if(result)
            message (STATUS "${ColorRed}Cannot configure ${LIBNAME}: ${result}${ColorEnd}")
            return()
        endif()

        execute_process( COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${LIBNAME}Download
        )

        if(result)
            message (STATUS "${ColorRed}Cannot build ${LIBNAME}: ${result}${ColorEnd}")
            return()
        endif()

        # If the header still does not exists, something went wrong.
        IF( NOT EXISTS ${LIBPATH}/${LIBNAME}/${TESTFILE} )
            message (STATUS "${ColorRed}Downloading ${LIBNAME} failed for unknown reasons${ColorEnd}")
            return()
        ENDIF()

        message (STATUS "${ColorBlue}Finished downloading ${LIBNAME}${ColorEnd}")
    ENDIF()

endfunction()

# This is a shortcut function for the above, which uses a github commit url to get the dependency.
function( DOWNLOAD_GITHUB_DEPENDENCY LIBPATH LIBNAME TESTFILE REPOSITORY COMMITHASH )

    set( LIBURL "https://github.com/${REPOSITORY}/archive/${COMMITHASH}.zip" )
    DOWNLOAD_DEPENDENCY( ${LIBPATH} ${LIBNAME} ${TESTFILE} ${LIBURL} )

endfunction()
