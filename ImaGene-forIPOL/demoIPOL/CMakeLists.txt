# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

set( DEMO_EXEC_FILES 
     meaningfulScaleEstim 
     minimalCodeNoiseDetection)	


set( DEMO_EXEC_FILES_EXT )
foreach( X ${DEMO_EXEC_FILES} )
	 set( DEMO_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${DEMO_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${DEMO_EXEC_NAME} ${LIBIMAGENE_NAME} )
	 list(APPEND BIN_EXEC_FILES_EXT ${DEMO_EXEC_NAME})
endforeach( X )

install( TARGETS ${DEMO_EXEC_FILES_EXT}
	 RUNTIME DESTINATION demoIPOL )

install( PROGRAMS ${SH_FILES} 
	 DESTINATION demoIPOL )

