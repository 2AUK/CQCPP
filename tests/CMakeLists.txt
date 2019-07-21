include_directories(../include)
add_executable(test_overlap test_integrals.cpp ${SOURCES})

add_test(overlap_test test_overlap)
