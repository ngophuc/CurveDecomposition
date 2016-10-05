# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

set( TESTPAMI_EXEC_FILES 
	displayNoise
	displayBoxes
	displayFlatCurved
	composeImage
	compareShape
	compareCntDistance
	ajusteImages	
	affTangentes)



set( TESTPAMI_EXEC_FILES_EXT )
foreach( X ${TESTPAMI_EXEC_FILES} )
	 set( TESTPAMI_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${TESTPAMI_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${TESTPAMI_EXEC_NAME} ${LIBIMAGENE_NAME} )
	 list(APPEND TEST_EXEC_FILES_EXT ${TESTPAMI_EXEC_NAME})
endforeach( X )









