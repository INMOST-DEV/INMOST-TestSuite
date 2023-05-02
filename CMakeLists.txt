cmake_minimum_required(VERSION 3.2)
project(TestSuite)

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})
#message(${INMOST_INCLUDE_DIRS})

#Navier-Stokes tests

add_executable(test_ns_linear ns_linear.cpp)
add_executable(test_ns_cavity ns_cavity.cpp)
add_executable(test_ns_pousielle ns_pousielle.cpp)
add_executable(test_ns_backward_step ns_backward_step.cpp)
add_executable(test_ns_turek ns_turek.cpp)
add_executable(test_ns_ethierstienman ns_ethiersteinman.cpp)
add_executable(test_ns_crevice ns_crevice.cpp)
add_executable(test_ns_ventrapp ns_ventrapp.cpp)

list(APPEND EXECS 
	test_ns_linear 
	test_ns_cavity 
	test_ns_pousielle
	test_ns_backward_step
	test_ns_turek
	test_ns_ethierstienman
	test_ns_crevice
	test_ns_ventrapp)
	
#Anisotropic diffusion tests

add_executable(test_diff_single_well single_well.cpp)
add_executable(test_diff_circular_perm circular_perm.cpp)
add_executable(test_diff_zigzag_perm zigzag_perm.cpp)
add_executable(test_diff_randmap_perm randmap_perm.cpp)
add_executable(test_diff_lognorm_perm lognorm_perm.cpp)
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

list(APPEND EXECS 
	test_diff_single_well 
	test_diff_circular_perm 
	test_diff_zigzag_perm 
	test_diff_randmap_perm 
	test_diff_lognorm_perm 
	test_diff_channels_perm
	test_diff_locking_perm
	test_diff_fvca5_test1_1
	test_diff_fvca5_test1_2
	test_diff_fvca5_test3
	test_diff_fvca5_test4
	test_diff_fvca5_test5
	test_diff_fvca5_test6
	test_diff_fvca5_test7
	test_diff_fvca5_test8
	test_diff_fvca5_test9
	test_diff_wugao_test2
	test_diff_wugao_test3
	test_diff_wugao_test4
	test_diff_two_wells_3d
	test_diff_edwards_3dtest4
	test_diff_discontinuous_4zones
	test_diff_linear_sol)

# Advection-diffusion tests

add_executable(test_advdiff_test00 adv00.cpp)
add_executable(test_advdiff_test01 adv01.cpp)
add_executable(test_advdiff_test02 adv02.cpp)
add_executable(test_advdiff_test03 adv03.cpp)
add_executable(test_advdiff_fvca5_test7 fvca5_test7_with_advection.cpp)
add_executable(test_advdiff_wugao_test3 wugao_test3_unsymmetric_with_advection.cpp)

list(APPEND EXECS 
	test_advdiff_test00 
	test_advdiff_test01 
	test_advdiff_test02
	test_advdiff_test03
	test_advdiff_fvca5_test7
	test_advdiff_wugao_test3)

# Tests for linear elasticity

add_executable(test_elastic_analytic1 elastic_analytic1.cpp)
add_executable(test_elastic_analytic2 elastic_analytic2.cpp)
add_executable(test_elastic_analytic3 elastic_analytic3.cpp)
add_executable(test_elastic_lshaped elastic_lshaped.cpp)
add_executable(test_elastic_membrane1 elastic_membrane1.cpp)
add_executable(test_elastic_membrane2 elastic_membrane2.cpp)
add_executable(test_elastic_infsphere elastic_infsphere.cpp)
add_executable(test_elastic_perforated elastic_perforated.cpp)
add_executable(test_elastic_beam_shear elastic_beam_shear.cpp)
add_executable(test_elastic_beam_torsion elastic_beam_torsion.cpp)
add_executable(test_elastic_contact1 elastic_contact1.cpp)
add_executable(test_elastic_contact2 elastic_contact2.cpp)

list(APPEND EXECS 
	test_elastic_analytic1 
	test_elastic_analytic2 
	test_elastic_analytic3
	test_elastic_lshaped
	test_elastic_membrane1
	test_elastic_membrane2
	test_elastic_infsphere
	test_elastic_perforated
	test_elastic_beam_shear
	test_elastic_beam_torsion
	test_elastic_contact1
	test_elastic_contact2)

#manipulation of grid and tensor

add_executable(test_transform_grid transform_grid.cpp)

list(APPEND EXECS 
	test_transform_grid)

foreach(EXEC ${EXECS})
	target_link_libraries(${EXEC} ${INMOST_LIBRARIES})
	install(TARGETS ${EXEC} EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)
endforeach(EXEC ${EXECS})



#poromechanics
add_subdirectory(poromech_test)

install(FILES readme.txt DESTINATION bin/TestSuite)

