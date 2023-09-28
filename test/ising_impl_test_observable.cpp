#include <vector>
#include <iostream>
#include <cassert>

#include "../include/ising_impl.hpp"
#include "../extern/mc_lib/mc_lib/_observable/observable.h"


void test()
{
    const uint32_t L = 16;
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


    IsingSimImpl isim(
        nei,
        nei_start,
        1.0,
        1.0,
        123);


    const uint32_t burn_in = 50000;
    const uint32_t n_iterations = 100000;

    struct agg_t {
        mc_stats::ScalarObservable<double> mag;
        mc_stats::ScalarObservable<double> mag_a;
        void operator()(int_fast32_t m, double e) {
            mag << m;
            mag_a << std::abs(m);
        };
    } agg;

    isim.iterate(n_iterations, agg, burn_in, 1);

    double mag_mean = agg.mag.mean() / N;
    double mag_a_mean = agg.mag_a.mean() / N;

    std::cout << "first measure\n";
    std::cout << "mag: " << mag_mean << std::endl;
    std::cout << "mag_abs: " << mag_a_mean << std::endl;
    assert(mag_mean < 0.01);
    assert(mag_a_mean > 0.5);


}

int main() {
    test();
    return 0;
}