# CMake generated Testfile for 
# Source directory: /home2/echar/imagene/tests
# Build directory: /home2/echar/imagene/tests
# 
# This file replicates the SUBDIRS() and ADD_TEST() commands from the source
# tree CMakeLists.txt file, skipping any SUBDIRS() or ADD_TEST() commands
# that are excluded by CMake control structures, i.e. IF() commands.
ADD_TEST(test_AffineFunction "test_AffineFunction")
ADD_TEST(test_HashTable "test_HashTable")
ADD_TEST(test_rSet "test_rSet")
ADD_TEST(test_DLine "test_DLine")
ADD_TEST(test_Vector "test_Vector")
ADD_TEST(test_Freeman "test_Freeman")
ADD_TEST(test_MLP "test_MLP")
ADD_TEST(test_Measure "test_Measure")
ADD_TEST(test_Math "test_Math" "-slr")
ADD_TEST(test_Arithmetic "test_Arithmetic" "-exp_stats_partial_quotients")
ADD_TEST(K2Space-slinel1 "test_K2Space" "-input" "/home2/echar/imagene/tests/chain3.fc" "-trace" "0" "-timing" "0" "-test_slinel")
ADD_TEST(K2Space-slinel2 "test_K2Space" "-input" "/home2/echar/imagene/tests/chain.fc" "-trace" "0" "-timing" "0" "-test_slinel")
SUBDIRS(TestPAMI)
SUBDIRS(TestCompNoiseDetect)
