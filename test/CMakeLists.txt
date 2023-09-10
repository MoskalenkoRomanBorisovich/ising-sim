add_executable(ising_impl_test ising_impl_test_main.cpp)

set(TEST1_FILES ising_impl_test.cpp ising_impl_test.hpp)
target_sources(ising_impl_test PRIVATE ${TEST1_FILES})

# I'm using C++17 in the test
target_compile_features(ising_impl_test PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(ising_impl_test ising_lib)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME ising_impl_test COMMAND ising_impl_test) # Command can be a target


# run test on build
add_custom_command(
     TARGET ising_impl_test
     COMMENT "Run tests"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> -R "^ising_impl_test" --output-on-failures
)