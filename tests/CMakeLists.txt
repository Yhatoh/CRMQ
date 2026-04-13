add_executable(test_huff test_huff.cpp)
target_link_libraries(test_huff PRIVATE cct_lib libsais)

add_executable(test_cct test_cct.cpp)
target_link_libraries(test_cct PRIVATE cct_lib libsais)

add_executable(test_rmq test_rmq.cpp)
target_link_libraries(test_rmq PRIVATE cct_lib libsais)

add_executable(unit_test_cct unit_test_cct.cpp)
target_link_libraries(unit_test_cct PRIVATE gtest gtest_main cct_lib libsais)

add_executable(unit_test_rmq unit_test_rmq.cpp)
target_link_libraries(unit_test_rmq PRIVATE gtest gtest_main cct_lib libsais)
