# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)


INCLUDE(CTest)

add_subdirectory (TestPAMI)
add_subdirectory (TestCompNoiseDetect)

set( TEST_EXEC_FILES 
        test_AffineFunction
 	test_Analog 
 	test_Arithmetic
 	test_HashTable 
 	test_TangentialCover 
 	test_nD 
 	test_rSet 
 	test_DSegment 
 	test_DLine 
 	test_Curvature 
 	test_Vector 
 	test_Freeman 
 	test_CubicalComplex 
 	test_Reconstruction
	test_K2Space 
	test_Resolutions
	test_MSOnContours
	test_Multiscale 
	test_Multiscale_segment 
	test_ImageScaleAnalysis
	test_OrderedAlphabet
	test_MLP 
	test_DMLP 
	test_DMLP_flip
	test_GridCurve
	test_MultNormales
	test_Profiles
	test_Measure
#	test_Tangent
	test_ErrorMeasure
	test_Extract3DCC
	test_BlurredSegmentTgtCover
	test_SBFraction
)


set( TEST_EXEC_FILES_EXT )
foreach( X ${TEST_EXEC_FILES} )
	 set( TEST_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${TEST_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${TEST_EXEC_NAME} ${LIBIMAGENE_NAME} )
	 list(APPEND TEST_EXEC_FILES_EXT ${TEST_EXEC_NAME})
endforeach( X )

set(TEST_EXEC_SIMPLETEST
     test_AffineFunction
		 test_HashTable
	   test_rSet
		 test_DLine
     test_Vector
     test_Freeman
     test_MLP
		 test_Measure
)
	
foreach( X ${TEST_EXEC_SIMPLETEST} )
	 set( TEST_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_test(${TEST_EXEC_NAME} ${TEST_EXEC_NAME})
endforeach( X )


#------------SPECIFIC TESTs----------
add_test(test_Math${SUFFIXBIN} test_Math${SUFFIXBIN} -slr)
add_test(test_Arithmetic${SUFFIXBIN} test_Arithmetic${SUFFIXBIN} -exp_stats_partial_quotients)


add_test( K2Space-slinel1 test_K2Space -input ${PROJECT_SOURCE_DIR}/tests/chain3.fc
	   -trace ${TRACE_LEVEL}
	    -timing ${TIMING_LEVEL}
	  -test_slinel )
add_test( K2Space-slinel2 test_K2Space -input ${PROJECT_SOURCE_DIR}/tests/chain.fc
	   -trace ${TRACE_LEVEL}
	    -timing ${TIMING_LEVEL}
	  -test_slinel )


