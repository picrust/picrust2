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

# ------------------------------------------------------------------------------
#   Info
# ------------------------------------------------------------------------------

#     This script is used for downloading a dependency if it is not found.
#     It is used as a build-time execution script by CMake.
#
#     We use @-substitution of the configure_file command of CMake
#     for the variables that determine which lib to download.
#     Those have to be set by the CMake script that uses this file:
#
#     DEPENDENCY_PATH: Path to the library dir, where the dependency is to be stored.
#     DEPENDENCY_NAME: Name of the dependency, determining the directory to download to.
#     DEPENDENCY_URL: The url of the file to download/extract.
#
#     In total, the download will thus end up in DEPENDENCY_PATH/DEPENDENCY_NAME

# ------------------------------------------------------------------------------
#   Download
# ------------------------------------------------------------------------------

# This min requirement is less than what we expect in the main CMakeList file,
# so we should be good. We state it here for re-use of this script.
cmake_minimum_required( VERSION 2.8.2 )

project( @DEPENDENCY_NAME@-download NONE )
include(ExternalProject)

# The download progress is ugly and not needed. Since CMake 3.1, we can disable it.
IF( ${CMAKE_VERSION} VERSION_GREATER 3.1 )
    SET( CMAKE_DOWNLOAD_PROGRESS "DOWNLOAD_NO_PROGRESS 1" )
ENDIF()

# message (STATUS "Downloading @DEPENDENCY_NAME@ from @DEPENDENCY_URL@")

# Download a fixed commit instead of the current master, so that we know that it works for us.
ExternalProject_Add( @DEPENDENCY_NAME@
    URL @DEPENDENCY_URL@
    SOURCE_DIR        "@DEPENDENCY_PATH@/@DEPENDENCY_NAME@"
    BINARY_DIR        "@DEPENDENCY_PATH@/@DEPENDENCY_NAME@"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
)
