# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
    endif(COMMAND cmake_policy)


INCLUDE(CTest)

set( EXAMPLES_EXEC_FILES 
     example-1-1
     example-1-2 
     example-1-3 
     example-1-3b 
     example-1-4
     example-1-5
      )

# set( EXAMPLES_EXEC_FILES_EXT )
foreach( X ${EXAMPLES_EXEC_FILES} )
	 set( EXAMPLES_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${EXAMPLES_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${EXAMPLES_EXEC_NAME} ${LIBIMAGENE_NAME} )
endforeach( X )


foreach( X ${EXAMPLES_EXEC_FILES_WITH_GMP} )
	 set( EXAMPLES_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${EXAMPLES_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${EXAMPLES_EXEC_NAME} ${LIBIMAGENE_NAME} ${GMP_LIBRARY} ${GMPXX_LIBRARY} )
endforeach( X )

set( EXAMPLES_EXEC_SIMPLETEST 
     example-1-1
     example-1-2 
     example-1-3 
     example-1-3b 
     example-1-4
     example-1-5
      )

foreach( X ${EXAMPLES_EXEC_SIMPLETEST} )
	 add_test(${EXAMPLES_EXEC_NAME} ${X})
endforeach( X )


# install( TARGETS ${EXAMPLES_EXEC_FILES_EXT}
# 	 RUNTIME DESTINATION bin )
# foreach( X ${EXAMPLES_EXEC_FILES} )
# 	 set( EXAMPLES_EXEC_NAME "${X}${SUFFIXBIN}" )
# 	 add_executable( ${EXAMPLES_EXEC_NAME} ${X}.cxx )
# 	 target_link_libraries( ${EXAMPLES_EXEC_NAME} ${LIBIMAGENE_NAME} )
# 	 list(APPEND EXAMPLES_EXEC_FILES_EXT ${EXAMPLES_EXEC_NAME})
# endforeach( X )
# target_link_libraries( exempleCOBA gmp gmpxx )

	 

