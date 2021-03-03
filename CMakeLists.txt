cmake_minimum_required(VERSION 2.8)
project(TestSuite)

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})
message(${INMOST_INCLUDE_DIRS})

#Navier-Stokes tests

add_executable(test_ns_linear ns_linear.cpp)
add_executable(test_ns_cavity ns_cavity.cpp)
add_executable(test_ns_pousielle ns_pousielle.cpp)
add_executable(test_ns_backward_step ns_backward_step.cpp)
add_executable(test_ns_turek ns_turek.cpp)
add_executable(test_ns_crevice ns_crevice.cpp)

target_link_libraries(test_ns_linear ${INMOST_LIBRARIES})
target_link_libraries(test_ns_cavity ${INMOST_LIBRARIES})
target_link_libraries(test_ns_pousielle ${INMOST_LIBRARIES})
target_link_libraries(test_ns_turek ${INMOST_LIBRARIES})
target_link_libraries(test_ns_crevice ${INMOST_LIBRARIES})
target_link_libraries(test_ns_backward_step ${INMOST_LIBRARIES})

install(TARGETS test_ns_linear EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_ns_cavity EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_ns_pousielle EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_ns_turek EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_ns_crevice EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_ns_backward_step EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

#Anisotropic diffusion tests

add_executable(test_diff_single_well single_well.cpp)
add_executable(test_diff_circular_perm circular_perm.cpp)
add_executable(test_diff_channels_perm channels_perm.cpp)
add_executable(test_diff_locking_perm locking_perm.cpp)
add_executable(test_diff_fvca5_test1_1 fvca5_test1_1.cpp)
add_executable(test_diff_fvca5_test1_2 fvca5_test1_2.cpp)
add_executable(test_diff_fvca5_test3 fvca5_test3.cpp)
add_executable(test_diff_fvca5_test4 fvca5_test4.cpp)
add_executable(test_diff_fvca5_test5 fvca5_test5.cpp)
add_executable(test_diff_fvca5_test6 fvca5_test6.cpp)
add_executable(test_diff_fvca5_test7 fvca5_test7.cpp)
add_executable(test_diff_fvca5_test8 fvca5_test8.cpp)
add_executable(test_diff_fvca5_test9 fvca5_test9.cpp)
add_executable(test_diff_wugao_test2 wugao_test2.cpp)
add_executable(test_diff_wugao_test3 wugao_test3_unsymmetric.cpp)
add_executable(test_diff_wugao_test4 wugao_test4.cpp)
add_executable(test_diff_two_wells_3d two_wells_3d.cpp)
add_executable(test_diff_edwards_3dtest4 edwards_3dtest4.cpp)
add_executable(test_diff_discontinuous_4zones discontinuous_4zones.cpp)
add_executable(test_diff_linear_sol linear_sol.cpp)

target_link_libraries(test_diff_single_well ${INMOST_LIBRARIES})
target_link_libraries(test_diff_circular_perm ${INMOST_LIBRARIES})
target_link_libraries(test_diff_channels_perm ${INMOST_LIBRARIES})
target_link_libraries(test_diff_locking_perm ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test1_1 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test1_2 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test3 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test4 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test5 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test6 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test7 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test8 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_fvca5_test9 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_wugao_test2 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_wugao_test3 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_wugao_test4 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_two_wells_3d ${INMOST_LIBRARIES})
target_link_libraries(test_diff_edwards_3dtest4 ${INMOST_LIBRARIES})
target_link_libraries(test_diff_discontinuous_4zones ${INMOST_LIBRARIES})
target_link_libraries(test_diff_linear_sol ${INMOST_LIBRARIES})


install(TARGETS test_diff_single_well EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_circular_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_channels_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_locking_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test1_1 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test1_2 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test3 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test5 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test6 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test7 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test8 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_fvca5_test9 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_wugao_test2 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_wugao_test3 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_wugao_test4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_two_wells_3d EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_edwards_3dtest4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_discontinuous_4zones EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_diff_linear_sol EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

# Advection-diffusion tests

add_executable(test_advdiff_test00 adv00.cpp)
add_executable(test_advdiff_test01 adv01.cpp)
add_executable(test_advdiff_test02 adv02.cpp)
add_executable(test_advdiff_test03 adv03.cpp)
add_executable(test_advdiff_fvca5_test7 fvca5_test7_with_advection.cpp)
add_executable(test_advdiff_wugao_test3 wugao_test3_unsymmetric_with_advection.cpp)

target_link_libraries(test_advdiff_test00 ${INMOST_LIBRARIES})
target_link_libraries(test_advdiff_test01 ${INMOST_LIBRARIES})
target_link_libraries(test_advdiff_test02 ${INMOST_LIBRARIES})
target_link_libraries(test_advdiff_test03 ${INMOST_LIBRARIES})
target_link_libraries(test_advdiff_fvca5_test7 ${INMOST_LIBRARIES})
target_link_libraries(test_advdiff_wugao_test3 ${INMOST_LIBRARIES})

install(TARGETS test_advdiff_test00 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_advdiff_test01 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_advdiff_test02 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_advdiff_test03 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_advdiff_fvca5_test7 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
install(TARGETS test_advdiff_wugao_test3 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


#manipulation of grid and tensor

add_executable(test_transform_grid transform_grid.cpp)

target_link_libraries(test_transform_grid ${INMOST_LIBRARIES})
install(TARGETS test_transform_grid EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)





install(FILES readme.txt DESTINATION bin/TestSuite)

