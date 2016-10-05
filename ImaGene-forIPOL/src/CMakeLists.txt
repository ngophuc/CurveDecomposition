# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)


set( ImaGene_SRC 
	base/Arguments.cxx 
	base/HashTable.cxx 
	base/Matrix.cxx 
	base/OrderedAlphabet.cxx
	base/Proxy.cxx 
	base/Signal.cxx 
	base/StandardArguments.cxx 
	base/StandardHashTable.cxx 
	base/StaticHashTable.cxx 
	base/RnMaps.cxx 
	base/Shapes.cxx 
	base/SmartSet.cxx 
	base/UndoableList.cxx 
	base/Vector.cxx 
	base/Vector2i.cxx 
	base/VectorUtils.cxx 
	base/WordUtils.cxx 
	timetools/Alarm.cxx 
	timetools/Clock.cxx 
        dgeometry2d/C4CELength.cxx 
        dgeometry2d/C4CGeometry.cxx 
        dgeometry2d/C4CIterator.cxx 
        dgeometry2d/C4CIteratorOnFreemanChain.cxx 
        dgeometry2d/C4CSegment.cxx 
        dgeometry2d/C4CSegmentPencil.cxx 
        dgeometry2d/C4CTangentialCover.cxx 
        dgeometry2d/CurveCode.cxx 
        dgeometry2d/DLine.cxx 
        dgeometry2d/DMLPContour.cxx 
        dgeometry2d/EuclideanGeometry.cxx 
	dgeometry2d/FreemanChain.cxx 
	dgeometry2d/FreemanChainTransform.cxx 
        dgeometry2d/FreemanFrame2D.cxx 
        dgeometry2d/GridCurve.cxx 
	dgeometry2d/MultiscaleFreemanChain.cxx 
	dgeometry2d/OctantChange.cxx 
	dgeometry2d/BlurredSegmentTgtCover.cxx 
	digitalnD/AreaOnSurfaceComputerByEuler.cxx  
        digitalnD/BelAdjacency.cxx 
        digitalnD/Bitset1.cxx 
        digitalnD/BoxEmbedder.cxx 
        digitalnD/C4CIteratorOnBdry.cxx 
        digitalnD/C4CIteratorOnExplicitDigitalSurface.cxx 
        digitalnD/C4CIteratorOnFreemanChainSurface.cxx 
	digitalnD/CentroidEmbedder.cxx 
	digitalnD/CubicalCell.cxx 
	digitalnD/CubicalComplex.cxx 
	digitalnD/CubicalSpace.cxx 
	digitalnD/DigitalSurface.cxx 
	digitalnD/DigitalSurfaceGeometry.cxx 
	digitalnD/DigitalSurfaceUmbrella.cxx 
	digitalnD/DigitalSurfaceUmbrellaComputer.cxx 
	digitalnD/Embedder.cxx 
	digitalnD/ExplicitDigitalSurface.cxx 
	digitalnD/ExplicitDigitalSurfaceTracker.cxx 
	digitalnD/Frame2D.cxx 
        digitalnD/GeometryComputerByContour4.cxx 
	digitalnD/GridEmbedder.cxx 
        digitalnD/K2Space.cxx 
        digitalnD/K3Isosurface.cxx 
        digitalnD/K3Isosurface_k6_l18.cxx 
        digitalnD/K3Shapes.cxx 
        digitalnD/KnCharSet.cxx 
        digitalnD/KnRCellSet.cxx 
        digitalnD/KnRCellVector.cxx 
        digitalnD/KnShapes.cxx 
        digitalnD/KnSpace.cxx 
        digitalnD/KnSpaceCoordScanner.cxx 
        digitalnD/KnSpaceScanner.cxx 
        digitalnD/KnTools.cxx 
        digitalnD/KnTypes.cxx 
        digitalnD/KnUtils.cxx 
        digitalnD/MeanCurvatureOnSurfaceComputerByDG.cxx 
        digitalnD/ObjectBoundary.cxx 
        digitalnD/ObjectBoundaryTracker.cxx 
        digitalnD/Sanitizer.cxx 
	digitalnD/SmoothingEmbedder.cxx 
        digitalnD/SurfelNeighborhood.cxx 
	digitalnD/NormalOnSurfaceComputerByUmbrellas.cxx	 
	helper/C4CTangentialCoverGeometry.cxx 
	helper/ChangeQuadrant.cxx 
	helper/CharacteristicPolygon.cxx 
	helper/ContourHelper.cxx 
	helper/CurveVariationsHelper.cxx 
	helper/DigitalElevationMap.cxx 
	helper/DigitalElevationMapStats.cxx 
	helper/DrawingXFIG.cxx
	helper/GlobalC4CGeometry.cxx 
	helper/ImageMultiscaleAnalysis.cxx 
	helper/LinearMinimizer.cxx 
	helper/MultiscaleProfile.cxx 
	helper/ScaleProfile.cxx 
	helper/ShapeHelper.cxx 
	helper/VisualCurvature.cxx 
	helper/Color.cxx 
	helper/MeaningfulContinuation.cxx
	image/Pixel.cxx
	image/Image2D.cxx
	image/Image2DUtils.cxx
	image/PGMFilter.cxx
	mathutils/CFraction.cxx 
	mathutils/Functions.cxx 
	mathutils/G.cxx 
	mathutils/DG.cxx 
	mathutils/Lattice.cxx 
	mathutils/Line2D.cxx 
	mathutils/Mathutils.cxx 
	mathutils/Polynomial.cxx 
	mathutils/Polynomial3.cxx 
	mathutils/PWFAffineFunction.cxx 
	mathutils/SampledFunction.cxx 
	mathutils/SimpleLinearRegression.cxx 
	mathutils/Statistic.cxx 
	mathutils/Statistics.cxx 
	mathutils/SternBrocot.cxx 
	mathutils/MeasureOfStraightLines.cxx
	3d/C4CIteratorOnImplicitDigitalSurface.cxx 
	3d/ImplicitDigitalSurface.cxx 
	3d/ImplicitDigitalSurfaceTracker.cxx 
	3d/ImplicitExtension.cxx 
	3d/ImplicitSurface.cxx 
	3d/ImplicitSurfaceHelper.cxx 
	3d/ImplicitSurfaceProjector.cxx 
	3d/ImplicitSurfaces.cxx 
	3d/ImplicitSurfaces3D.cxx 
	3d/ImplicitSurfaceEmbedder.cxx 
        )


# Create a library called "ImaGene".
# The extension is already found. Any number of sources could be listed here.
add_library ( ${LIBIMAGENE_NAME} ${ImaGene_SRC}  )
target_link_libraries( ${LIBIMAGENE_NAME}
                       ${GMP_LIBRARY}
		       ${GMPXX_LIBRARY} )

install( TARGETS ${LIBIMAGENE_NAME} LIBRARY DESTINATION lib  )
