// #pragma once
#include <cstdint>
#include <optional>
#include <random>
#include <memory>
#include <cassert>

#include "ising_utils.hpp"


class IsingSimImpl
{
public:
    // typedefs
    using spin_t = int_fast8_t;
    using site_t = uint_fast32_t;
    typedef std::mt19937_64 Generator_t;

    template<typename S, typename P>
    IsingSimImpl(
        const std::vector<S>& neighbors,
        const std::vector<P>& nei_start,
        const double beta,
        const double J,
        const Generator_t::result_type seed = 12345,
        const std::vector<int>* cur_spins = nullptr);

    void iteration();

    // getters
    inline const std::vector<int8_t>& get_spins() const { return m_cur_spins; };
    inline const int_fast32_t get_mag() { return m_cur_mag; };
    inline const double get_ene() { return -m_J * m_cur_ene_spins; };
    inline const site_t get_N() { return m_N; };
    inline const double get_beta() { return m_beta; };
    inline const double get_J() { return m_J; };

private:
    const site_t m_N; // number of spins
    const site_t m_N_sqrt;
    const double m_beta;
    const double m_J;

    std::vector<site_t>  m_nei_start;
    std::vector<site_t> m_neighbors;

    void random_spins();

    const double m_accept_ratio; // cluster algorithm accept ration

    uint_fast32_t m_cur_step;

    std::vector<spin_t> m_cur_spins;

    Generator_t m_rng;
    std::uniform_real_distribution<double> m_dis;

    inline double m_rand() { return m_dis(m_rng); };

    // current model stats
    int_fast32_t m_cur_mag;
    int_fast32_t m_cur_ene_spins;
    void calc_cur_mag();
    void calc_cur_ene_spins(); // sum of spin[i] * spin[j]

};


// templates definition
template<typename S, typename P>
IsingSimImpl::IsingSimImpl(
    const std::vector<S>& neighbors,
    const std::vector<P>& nei_start,
    const double beta,
    const double J,
    const Generator_t::result_type seed,
    const std::vector<int>* cur_spins) :
    m_N(nei_start.size() - 1),
    m_N_sqrt(std::sqrt(m_N)),
    m_beta(beta),
    m_J(J),
    m_neighbors(neighbors.begin(), neighbors.end()),
    m_nei_start(nei_start.begin(), nei_start.end()),
    m_cur_step(0),
    // m_dis(0., 1.),
    m_accept_ratio(1.0 - std::exp(-2.0 * m_beta))
{
    assert(m_nei_start.size() == m_N + 1);
    assert(m_neighbors.size() == m_nei_start.back());
    sort_neighbors(m_nei_start, m_neighbors);
    m_rng.seed(seed);
    m_dis.reset();
    if (cur_spins) {
        assert(cur_spins->size() == m_N);
        m_cur_spins.clear();
        std::copy(cur_spins->begin(), cur_spins->end(), std::back_inserter(m_cur_spins));
    }
    else
        random_spins();

    calc_cur_ene_spins();
    calc_cur_mag();
}