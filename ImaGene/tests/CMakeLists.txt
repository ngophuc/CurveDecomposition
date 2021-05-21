# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)


INCLUDE(CTest)

add_subdirectory (TestCompNoiseDetect)

