project(POROMECH_TEST)

add_library(POROMECH_TEST 
		abstract_test.cpp abstract_test.h
		test1.cpp test1.h
		test2.cpp test2.h
		test3.cpp test3.h
		test4.cpp test4.h
		test5.cpp test5.h
		test6.cpp test6.h
		test7.cpp test7.h
        test8.cpp test8.h
		test10.cpp test10.h
		test11.cpp test11.h
		test13.cpp test13.h)

set_property(TARGET POROMECH_TEST PROPERTY PUBLIC_HEADER abstract_test.h)

target_include_directories(POROMECH_TEST
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

#target_include_directories(POROMECH_TEST INTERFACE .)

#install(FILES abstract_test.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS POROMECH_TEST EXPORT POROMECH_TEST-Targets 
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

install(EXPORT POROMECH_TEST-Targets
        FILE POROMECH_TEST-Targets.cmake
#		NAMESPACE POROMECH_TEST::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/POROMECH_TEST)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/POROMECH_TESTConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/POROMECH_TEST)

set(version 0.1.1)
set_property(TARGET POROMECH_TEST PROPERTY VERSION ${version})
set_property(TARGET POROMECH_TEST PROPERTY SOVERSION 3)
set_property(TARGET POROMECH_TEST PROPERTY INTERFACE_POROMECH_TEST_MAJOR_VERSION 1)
set_property(TARGET POROMECH_TEST APPEND PROPERTY COMPATIBLE_INTERFACE_STRING POROMECH_TEST_MAJOR_VERSION)


# generate the version file for the config file
write_basic_package_version_file( "${CMAKE_CURRENT_BINARY_DIR}/POROMECH_TESTConfigVersion.cmake"
  VERSION "${version}" COMPATIBILITY AnyNewerVersion)

install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/POROMECH_TESTConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/POROMECH_TESTConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/POROMECH_TEST)

export(EXPORT POROMECH_TEST-Targets 
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/POROMECH_TEST-Targets.cmake"
#       NAMESPACE POROMECH_TEST::
)

export(TARGETS POROMECH_TEST FILE POROMECH_TEST-Targets.cmake)
export(PACKAGE POROMECH_TEST)