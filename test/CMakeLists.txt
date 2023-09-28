add_executable(ising_impl_test ising_impl_test_main.cpp)
set(TEST1_FILES ising_impl_test.cpp)
target_sources(ising_impl_test PRIVATE ${TEST1_FILES})
target_link_libraries(ising_impl_test ising_lib)
add_test(NAME ising_impl_test COMMAND ising_impl_test)

add_executable(ising_impl_test_observable ising_impl_test_observable.cpp)
target_link_libraries(ising_impl_test_observable ising_lib)
add_test(NAME ising_impl_test_observable COMMAND ising_impl_test_observable)

add_executable(ising_impl_benchmark ising_impl_benchmark.cpp)
target_link_libraries(ising_impl_benchmark ising_lib)

# run test on release builds
IF(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    message("release mode")
    add_custom_command(
        TARGET ising_impl_test
        COMMENT "Run basic tests"
        POST_BUILD 
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> -R "^ising_impl_test" --output-on-failures
    )
ENDIF()