# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

set( BIN_EXEC_FILES 
     	addGaussianNoisePGM
     	addGaussianNoisePGM3d
	christoffel
	curvature_cc 
	curvature_conv 
	curvature_gmc 
	curvature_visualc 
	error_analysis 
	freeman2fig 
	freeman2freeman 
	freeman2pgm 
	freeman2polygon 
	freeman2sdp 
	gencontour 
	integrate_shape 
	length_tan
	pgm2freeman 
	rudin-noise-removal
	sample_curve
	scale_analysis
	sdp2pts
	polygonTransform 
	gen3dShapes
	pgm2pgm3d)	

set ( SH_FILES
    analyse-CC-GMC-integration.sh
    analyse-error-CC-GMC.sh
    analyse-GMC-integration-close.sh
    maximal-segments.sh
    analyse-GMC-integration.sh
    polygon-digital-visualc.sh
)

set( BIN_EXEC_FILES_EXT )
foreach( X ${BIN_EXEC_FILES} )
	 set( BIN_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${BIN_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${BIN_EXEC_NAME} ${LIBIMAGENE_NAME} )
	 list(APPEND BIN_EXEC_FILES_EXT ${BIN_EXEC_NAME})
endforeach( X )

install( TARGETS ${BIN_EXEC_FILES_EXT}
	 RUNTIME DESTINATION bin )

install( PROGRAMS ${SH_FILES} 
	 DESTINATION bin )

