add_executable(test_huff test_huff.cpp)
target_link_libraries(test_huff PRIVATE CRMQ_lib libsais)

add_executable(test_scst test_scst.cpp)
target_link_libraries(test_scst PRIVATE CRMQ_lib libsais)

add_executable(test_rmq test_rmq.cpp)
target_link_libraries(test_rmq PRIVATE CRMQ_lib libsais)

add_executable(unit_test_scst unit_test_scst.cpp)
target_link_libraries(unit_test_scst PRIVATE gtest gtest_main CRMQ_lib libsais)

add_executable(unit_test_rmq unit_test_rmq.cpp)
target_link_libraries(unit_test_rmq PRIVATE gtest gtest_main CRMQ_lib libsais)
