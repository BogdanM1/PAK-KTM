cmake_minimum_required (VERSION 2.6)
project (pakktm Fortran CXX C)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
message( STATUS "WIN32: " ${WIN32} )
message( STATUS "UNIX:  " ${UNIX} )
message( STATUS "HOME:  " $ENV{HOME} )
message( STATUS "PROJECT_BINARY_DIR:  " ${PROJECT_BINARY_DIR} )
message( STATUS "PROJECT_SOURCE_DIR:  " ${PROJECT_SOURCE_DIR} )
message( STATUS "PROJECT_NAME:        " ${PROJECT_NAME} )
message( STATUS "CMAKE_CONFIGURATION_TYPES: " ${CMAKE_CONFIGURATION_TYPES} )
message( STATUS "CMAKE_MODULE_PATH: "         ${CMAKE_MODULE_PATH} )


# portable linux exe
set(CMAKE_INSTALL_RPATH "\$ORIGIN/libs")
set(CMAKE_BUILD_WITH_INSTALL_RPATH on)

file(GLOB_RECURSE PAK_SRC_FILES "src/pak/*.F90" "src/pak/*.f90" "src/pak/*.for" "src/pak/*.FOR" 
                                "src/libseq/*.f" "src/libseq/*.c"
                                "src/sparse_assembler/*.c")							
add_executable(pakktm ${PAK_SRC_FILES})

set_target_properties(pakktm PROPERTIES LINKER_LANGUAGE Fortran) 
set(CMAKE_Fortran_FLAGS "-cpp -g -O0 ${CMAKE_Fortran_FLAGS}") 

include_directories("src/MUMPS" "src/libseq")
TARGET_LINK_LIBRARIES(pakktm "${CMAKE_CURRENT_BINARY_DIR}/lib/libmumps_common.so"      
                           "${CMAKE_CURRENT_BINARY_DIR}/lib/libdmumps.so" 
                           "${CMAKE_CURRENT_BINARY_DIR}/lib/libpord.so")

message( STATUS "CMAKE_Fortran_FLAGS: "      ${CMAKE_Fortran_FLAGS})			
message( STATUS "CMAKE_Fortran_COMPILER: "   ${CMAKE_Fortran_COMPILER})



