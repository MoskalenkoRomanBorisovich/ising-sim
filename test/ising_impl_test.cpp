#include <vector>
#include <iostream>
#include <cassert>

#include "ising_impl_test.hpp"

#include "../include/ising_impl.hpp"


/*
basic test
*/
void test_1()
{
    const std::vector<int> nei = { 1, 2, 0, 2, 0, 1 };
    const std::vector<int> nei_start = { 0, 2, 4, 6 };


    IsingSimImpl isim(
        nei,
        nei_start,
        1.0,
        1.0,
        123);

    const int m1 = isim.get_mag();
    std::cout << m1 << "\n";
    isim.iteration();

    const int m2 = isim.get_mag();
    std::cout << m2 << "\n";

    assert(m1 != m2);
}

/*
basic 1D test
*/
void test_2()
{
    const uint32_t N = 250;

    std::vector<IsingSimImpl::site_t> nei = { 1, N - 1 };
    std::vector<IsingSimImpl::site_t> nei_start = { 0, 2 };
    for (uint32_t i = 1; i < N - 1; ++i) {
        nei.push_back(i - 1);
        nei.push_back(i + 1);
        nei_start.push_back(i * 2 + 2);
    }
    nei.push_back(N - 2);
    nei.push_back(0);
    nei_start.push_back(N * 2);

    IsingSimImpl isim(
        nei,
        nei_start,
        1.0,
        1.0,
        123);

    const uint32_t max_it = 100000;
    for (uint32_t i = 0; i < max_it; ++i) { // burn in
        isim.iteration();
    }
    double mag_mean = 0.0;
    for (uint32_t i = 0; i < max_it; ++i) {
        isim.iteration();
        mag_mean += std::abs(isim.get_mag()) / (double)N;
    }
    mag_mean /= max_it;
    std::cout << mag_mean << std::endl;
    assert(mag_mean < 0.5);
}