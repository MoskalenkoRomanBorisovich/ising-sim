#include <vector>
#include <iostream>
#include <cassert>

#include "../include/ising_impl.hpp"


void test_square32()
{
    std::cout << "Run test_square32\n";
    const uint32_t L = 32;
    const uint32_t N = L * L;


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


    IsingSimImpl isim;
    isim.init(
        nei,
        nei_start,
        1.0,
        1.0,
        123);


    const uint32_t burn_in = 100000;
    const uint32_t n_iterations = 500000;

    long long mag = 0;
    uint32_t mag_a = 0;
    const auto& agg = [&mag, &mag_a](int mag_, double ene_) {
        mag += mag_;
        mag_a += std::abs(mag_);
        return;
        };

    clock_t start = clock();
    isim.burn_in(burn_in);
    isim.iterate(n_iterations, agg, 1);
    clock_t end = clock();
    std::cout << "time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

    double mag_mean = mag / (double)(N * n_iterations);
    double mag_abs_mean = mag_a / (double)(N * n_iterations);
    std::cout << "first measure\n";
    std::cout << "mag: " << mag_mean << std::endl;
    std::cout << "mag_abs: " << mag_abs_mean << std::endl;
    assert(std::abs(mag_mean) < 0.01);
    assert(mag_abs_mean > 0.5);
}


void benchmark_aggregators()
{
    std::cout << "Run benchmark_aggregators\n";
    const uint32_t L = 32;
    const uint32_t N = L * L;


    std::vector<IsingSimImpl::site_t> nei;
    nei.reserve(N * 4);
    std::vector<IsingSimImpl::site_t> nei_start;
    nei_start.reserve(N + 1);
    for (uint32_t i = 0; i < L; ++i) {
        for (uint32_t j = 0; j < L; ++j) {
            const uint32_t idx = i * L + j;
            nei.push_back((N + idx - L) % N);
            nei.push_back((N + idx + L) % N);
            nei.push_back((N + idx - 1) % N);
            nei.push_back((N + idx + 1) % N);
            nei_start.push_back(idx * 4);
        }
    }
    nei_start.push_back(N * 4);

    const uint32_t burn_in = 0;
    const uint64_t n_iterations = 100000000L;
    const double beta = 0.1;
    { // lambda function
        IsingSimImpl isim;
        isim.init(
            nei,
            nei_start,
            beta,
            1.0,
            123);

        long long mag = 0;
        uint32_t mag_a = 0;
        const auto& agg = [&mag, &mag_a](int mag_, double ene_) {
            mag += mag_;
            mag_a += std::abs(mag_);
            };

        clock_t start = clock();
        isim.burn_in(burn_in);
        isim.iterate(n_iterations, agg, 1);
        clock_t end = clock();
        std::cout << "lambda time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

        double mag_mean = mag / (double)(N * n_iterations);
        double mag_abs_mean = mag_a / (double)(N * n_iterations);
        std::cout << "first measure\n";
        std::cout << "mag: " << mag_mean << std::endl;
        std::cout << "mag_abs: " << mag_abs_mean << std::endl;
        assert(std::abs(mag_mean) < 0.01);
        assert(mag_abs_mean > 0.5);
    }

    { // struct
        IsingSimImpl isim;
        isim.init(
            nei,
            nei_start,
            beta,
            1.0,
            123);

        struct my_agg {
            long long mag = 0;
            uint32_t mag_a = 0;
            inline void operator()(int mag_, double ene_) {
                mag += mag_;
                mag_a += std::abs(mag_);
            }
        } agg;

        clock_t start = clock();
        isim.burn_in(burn_in);
        isim.iterate(n_iterations, agg, 1);
        clock_t end = clock();
        std::cout << "struct time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

        double mag_mean = agg.mag / (double)(N * n_iterations);
        double mag_abs_mean = agg.mag_a / (double)(N * n_iterations);
        std::cout << "first measure\n";
        std::cout << "mag: " << mag_mean << std::endl;
        std::cout << "mag_abs: " << mag_abs_mean << std::endl;
        assert(std::abs(mag_mean) < 0.01);
        assert(mag_abs_mean > 0.5);
    }

    { // straight 
        IsingSimImpl isim;
        isim.init(
            nei,
            nei_start,
            beta,
            1.0,
            123);

        long long mag = 0;
        uint32_t mag_a = 0;

        clock_t start = clock();
        isim.burn_in(burn_in);
        for (uint_fast32_t i = 0; i < n_iterations; ++i) {
            isim.iterate();
            mag += isim.get_mag();
            mag_a += std::abs(isim.get_mag());
        }
        clock_t end = clock();
        std::cout << "straight time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

        double mag_mean = mag / (double)(N * n_iterations);
        double mag_abs_mean = mag_a / (double)(N * n_iterations);
        std::cout << "first measure\n";
        std::cout << "mag: " << mag_mean << std::endl;
        std::cout << "mag_abs: " << mag_abs_mean << std::endl;
        assert(std::abs(mag_mean) < 0.01);
        assert(mag_abs_mean > 0.5);
    }


}

int main() {
    benchmark_aggregators();
    test_square32();
    return 0;
}