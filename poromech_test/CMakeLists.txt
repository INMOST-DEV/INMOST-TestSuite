project(POROMECH_TEST)

add_library(POROMECH_TEST 
		abstract_test.cpp 
		test1.cpp 
		test2.cpp 
		test3.cpp 
		test4.cpp 
		test6.cpp  
		test7.cpp 
		test10.cpp 
		test11.cpp 
		test13.cpp)

set_property(TARGET POROMECH_TEST PROPERTY PUBLIC_HEADER abstract_test.h)

target_include_directories(POROMECH_TEST INTERFACE .)

install(TARGETS POROMECH_TEST EXPORT projects-targets 
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

