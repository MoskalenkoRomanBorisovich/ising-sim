#pragma once
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
    using spin_t = signed char; // spin value type
    using site_t = uint_fast32_t; // site index type
    using iter_t = uint_fast64_t; // iteration index type
    typedef std::mt19937_64 Generator_t;

    template<typename S, typename P, typename T = spin_t>
    void init(
        const std::vector<S>& neighbors,
        const std::vector<P>& nei_start,
        const double beta,
        const double J,
        const Generator_t::result_type seed = 12345,
        const std::vector<T>& cur_spins = std::vector<int_fast8_t>());

    void iterate();

    inline void iterate(iter_t n_steps)
    {
        for (iter_t i = 0; i < n_steps; ++i) {
            iterate();
        }
    }

    inline void burn_in(iter_t n_steps) // explicit alias for iterate
    {
        iterate(n_steps);
    }

    template<typename T>
    void iterate(iter_t n_measures, T& aggregator, iter_t measure_steps = 1);

    // getters
    inline const std::vector<int8_t>& get_spins() const { return m_cur_spins; };
    inline const int32_t get_mag() const { return m_cur_mag; };
    inline const double get_ene() const { return -m_J * m_cur_ene_spins; };
    inline const uint32_t get_N() const { return m_N; };
    inline const double get_beta() const { return m_beta; };
    inline const double get_J() const { return m_J; };

    void clear_buf() {
        m_cluster.clear();
        m_cluster.shrink_to_fit();
    }

private:
    site_t m_N; // number of spins
    double m_beta;
    double m_J;

    std::vector<site_t>  m_nei_start;
    std::vector<site_t> m_neighbors;

    void random_spins();

    double m_accept_ratio; // cluster algorithm accept ratio

    uint_fast32_t m_cur_step;

    std::vector<spin_t> m_cur_spins;

    Generator_t m_rng;
    std::uniform_real_distribution<double> m_dis;

    inline double m_rand() { return m_dis(m_rng); }

    std::vector<site_t> m_cluster; // buffer for cluster algorithm

    // current model stats
    int_fast32_t m_cur_mag;
    int_fast32_t m_cur_ene_spins;
    void calc_cur_mag();
    void calc_cur_ene_spins(); // sum of spin[i] * spin[j]

    template<spin_t spin_val>
    void recalculate_flipped_cluster_impl(const std::vector<site_t>& flipped_cluster);
    void recalculate_flipped_cluster(const std::vector<site_t>& flipped_cluster, const spin_t spin_val)
    {
        switch (spin_val)
        {
        case 1:
            recalculate_flipped_cluster_impl<1>(flipped_cluster);
            break;

        case -1:
            recalculate_flipped_cluster_impl<-1>(flipped_cluster);
            break;

        default:
            break;
        }
    }
};

// template functions
template<typename S, typename P, typename T = IsingSimImpl::spin_t>
void IsingSimImpl::init(
    const std::vector<S>& neighbors,
    const std::vector<P>& nei_start,
    const double beta,
    const double J,
    const Generator_t::result_type seed,
    const std::vector<T>& cur_spins)
{
    m_N = nei_start.size() - 1;
    m_beta = beta;
    m_J = J;
    m_neighbors = std::vector<site_t>(neighbors.begin(), neighbors.end());
    m_nei_start = std::vector<site_t>(nei_start.begin(), nei_start.end());
    m_cur_step = 0;
    m_dis = std::uniform_real_distribution<double>(0., 1.);
    m_accept_ratio = 1.0 - std::exp(-2.0 * m_beta);

    assert(m_nei_start.size() == m_N + 1);
    assert(m_neighbors.size() == m_nei_start.back());
    sort_neighbors(m_nei_start, m_neighbors);
    m_rng.seed(seed);
    m_dis.reset();
    if (!cur_spins.empty()) {
        assert(cur_spins.size() == m_N);
        m_cur_spins.clear();
        std::copy(cur_spins.begin(), cur_spins.end(), std::back_inserter(m_cur_spins));
    }
    else
        random_spins();

    calc_cur_ene_spins();
    calc_cur_mag();
}

template<typename T>
inline void IsingSimImpl::iterate(iter_t n_measures, T& aggregator, iter_t measure_steps)
{
    assert(0 < measure_steps);
    for (iter_t block = 0, block_end = n_measures; block < block_end; ++block) {
        iterate(measure_steps);
        aggregator(get_mag(), get_ene());
    }
}
