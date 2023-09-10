#include <vector>
#include <iostream>
#include <cassert>

#include "ising_impl_test.hpp"

#include "../include/ising_impl.hpp"



void test_1()
{

    std::cout << "Tests begin\n";
    const std::vector<IsingSimImpl::site_t> nei = { 1, 2, 0, 2, 0, 1 };
    const std::vector<IsingSimImpl::site_t> nei_start = { 0, 2, 4, 6 };


    IsingSimImpl isim(
        nei,
        nei_start,
        1.0,
        1.0,
        123,
        nullptr);

    const int m1 = isim.get_mag();
    std::cout << m1 << "\n";
    isim.iteration();

    const int m2 = isim.get_mag();
    std::cout << m2 << "\n";

    assert(m1 != m2);
}