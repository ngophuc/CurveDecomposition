# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

set( TESTCOMPNOISEDETECT_EXEC_FILES 
     displayNoiseBS  displayMeaningfulContours displayEpsilonMC testEpsilonMS )





set( TESTCOMPNOISEDETECT_EXEC_FILES_EXT )
foreach( X ${TESTCOMPNOISEDETECT_EXEC_FILES} )
	 set( TESTCOMPNOISEDETECT_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${TESTCOMPNOISEDETECT_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${TESTCOMPNOISEDETECT_EXEC_NAME} ${LIBIMAGENE_NAME} )
	 list(APPEND TEST_EXEC_FILES_EXT ${TESTCOMPNOISEDETECT_EXEC_NAME})
endforeach( X )









