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
    isim.iterate();

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
        isim.iterate();
    }
    double mag_mean = 0.0;
    for (uint32_t i = 0; i < max_it; ++i) {
        isim.iterate();
        mag_mean += std::abs(isim.get_mag()) / (double)N;
    }
    mag_mean /= max_it;
    std::cout << mag_mean << std::endl;
    assert(mag_mean < 0.5);
}

/*
    test iterate with aggrigator
*/
void test_3()
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

    const uint32_t burn_in = 50000;
    const uint32_t n_iterations = 100000;

    long long mag = 0;
    uint32_t mag_a = 0;
    const auto& agg = [&mag, &mag_a](int mag_, double ene_) {
        mag += mag_;
        mag_a += std::abs(mag_);
        return;
        };
    isim.iterate(n_iterations, agg, burn_in, 1);

    double mag_mean = mag / (double)(N * n_iterations);
    double mag_abs_mean = mag_a / (double)(N * n_iterations);
    std::cout << "first measure\n";
    std::cout << "mag: " << mag_mean << std::endl;
    std::cout << "mag_abs: " << mag_abs_mean << std::endl;
    assert(std::abs(mag_mean) < 0.01);
    assert(mag_abs_mean < 0.5);


    isim.iterate(n_iterations, agg);
    mag_mean = mag / (double)(N * n_iterations * 2);
    mag_abs_mean = mag_a / (double)(N * n_iterations * 2);
    std::cout << "second measure\n";
    std::cout << "mag: " << mag_mean << std::endl;
    std::cout << "mag_abs: " << mag_abs_mean << std::endl;
    assert(std::abs(mag_mean) < 0.01);
    assert(mag_abs_mean < 0.5);
}


/*
    test 2d grid
*/

void test_4()
{
    const uint32_t L = 10;
    const uint32_t N = 10 * 10;


    std::vector<IsingSimImpl::site_t> nei;
    nei.reserve(N * 4);
    std::vector<IsingSimImpl::site_t> nei_start;
    nei_start.reserve(N + 1);
    for (uint32_t i = 0; i < L; ++i) {
        for (uint32_t j = 0; j < L; ++j) {
            const uint32_t idx = i * L + j;
            nei.push_back((N + idx - L) % N);
            nei.push_back((N + idx - L) % N);
            nei.push_back((N + idx - 1) % N);
            nei.push_back((N + idx + 1) % N);
            nei_start.push_back(idx * 4);
        }
    }
    nei_start.push_back(N * 4);


    IsingSimImpl isim(
        nei,
        nei_start,
        1.0,
        1.0,
        123);


    const uint32_t burn_in = 50000;
    const uint32_t n_iterations = 100000;

    long long mag = 0;
    uint32_t mag_a = 0;
    const auto& agg = [&mag, &mag_a](int mag_, double ene_) {
        mag += mag_;
        mag_a += std::abs(mag_);
        return;
        };

    isim.iterate(n_iterations, agg, burn_in, 1);


    double mag_mean = mag / (double)(N * n_iterations);
    double mag_abs_mean = mag_a / (double)(N * n_iterations);
    std::cout << "first measure\n";
    std::cout << "mag: " << mag_mean << std::endl;
    std::cout << "mag_abs: " << mag_abs_mean << std::endl;
    assert(std::abs(mag_mean) < 0.01);
    assert(mag_abs_mean > 0.95);
}