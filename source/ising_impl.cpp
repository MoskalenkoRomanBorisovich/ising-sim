#include "../include/ising_impl.hpp"
#include "../include/ising_utils.hpp"
#include <unordered_set>


/*
    One step of cluster update
*/
void IsingSimImpl::iterate()
{
    std::vector<site_t> pocket;
    std::unordered_set<site_t> cluster;
    pocket.reserve(m_N_sqrt);
    pocket.push_back(m_rng() % m_N);
    cluster.insert(pocket.back());
    const spin_t spin_val = m_cur_spins[pocket.back()];
    const spin_t neg_spin_val = -spin_val;
    const spin_t neg_spin_val2 = neg_spin_val * 2;
    do {
        const site_t cur = pocket.back();
        pocket.pop_back();
        m_cur_spins[cur] = neg_spin_val;
        m_cur_mag += neg_spin_val2;
        for (site_t nei_id = m_nei_start[cur], nei_end = m_nei_start[cur + 1]; nei_id < nei_end; ++nei_id) {
            const site_t nei = m_neighbors[nei_id];
            if (spin_val == m_cur_spins[nei]) {
                m_cur_ene_spins -= 2;
                if (cluster.count(nei) == 0) {
                    if (m_rand() < m_accept_ratio) {
                        cluster.insert(nei);
                        pocket.push_back(nei);
                    }
                }
            }
            else {
                m_cur_ene_spins += 2;
            }
        }
    } while (!pocket.empty());
}

void IsingSimImpl::random_spins()
{
    m_cur_spins.resize(m_N);
    for (site_t i = 0; i < m_N; ++i) {
        m_cur_spins[i] = m_rand() < 0.5 ? -1 : 1;
    }
}

void IsingSimImpl::calc_cur_ene_spins()
{
    m_cur_ene_spins = 0; // reset energy
    for (site_t site = 0; site < m_N; ++site) {
        for (site_t nei = m_nei_start[site], nei_end = m_nei_start[site + 1]; nei < nei_end; ++nei) {
            const site_t site1 = m_neighbors[nei];
            if (site > site1) // neighbors are sorted, prevents double checking 
                m_cur_ene_spins += m_cur_spins[site] * m_cur_spins[site1];
            else
                break;
        }
    }
}

void IsingSimImpl::calc_cur_mag()
{
    m_cur_mag = std::accumulate(m_cur_spins.begin(), m_cur_spins.end(), 0);
}
