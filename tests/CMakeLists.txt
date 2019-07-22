include_directories(../include)
add_executable(test_overlap test_integrals.cpp ${SOURCES})

set(EXTERNAL_LIBS cblas lapack gsl)

target_link_libraries(test_overlap ${EXTERNAL_LIBS})

add_test(overlap_test test_overlap)
