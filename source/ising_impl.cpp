#include "../include/ising_impl.hpp"
#include "../include/ising_utils.hpp"
#include <unordered_set>
#include "ising_impl.hpp"


template<IsingSimImpl::spin_t spin_val>
void IsingSimImpl::recalculate_flipped_cluster_impl(const std::vector<site_t>& flipped_cluster)
{
    m_cur_mag -= 2 * spin_val * flipped_cluster.size();

    for (site_t cluster_pos = 0, cluster_end = flipped_cluster.size(); cluster_pos < cluster_end; ++cluster_pos) {
        const site_t cur = flipped_cluster[cluster_pos];
        for (site_t nei_id = m_nei_start[cur], nei_end = m_nei_start[cur + 1]; nei_id < nei_end; ++nei_id) {
            const site_t nei = m_neighbors[nei_id];
            switch (m_cur_spins[nei])
            {
            case spin_val:
                m_cur_ene_spins -= 2;
                break;

            case -spin_val:
                m_cur_ene_spins += 2;
                break;
            default:
                break;
            }
        }
    }
    for (site_t cluster_pos = 0, cluster_end = flipped_cluster.size(); cluster_pos < cluster_end; ++cluster_pos) {
        const site_t cur = flipped_cluster[cluster_pos];
        m_cur_spins[cur] = -spin_val;
    }
}
/*
    One step of cluster update
*/
void IsingSimImpl::iterate()
{
    m_cluster.clear();
    m_cluster.push_back(m_rng() % m_N);
    const spin_t spin_val = m_cur_spins[m_cluster[0]];
    m_cur_spins[m_cluster[0]] = 0;
    for (site_t cluster_pos = 0; cluster_pos < m_cluster.size(); ++cluster_pos) {
        const site_t cur = m_cluster[cluster_pos];
        for (site_t nei_id = m_nei_start[cur], nei_end = m_nei_start[cur + 1]; nei_id < nei_end; ++nei_id) {
            const site_t nei = m_neighbors[nei_id];
            if (spin_val == m_cur_spins[nei]) {
                if (m_rand() < m_accept_ratio) {
                    m_cluster.push_back(nei);
                    m_cur_spins[nei] = 0;
                }
            }
        }
    }
    recalculate_flipped_cluster(m_cluster, spin_val);
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
            //if (site > site1) // neighbors are sorted, prevents double checking 
            m_cur_ene_spins += m_cur_spins[site] * m_cur_spins[site1];
            //else
            //    break;
        }
    }
    m_cur_ene_spins /= 2;
}


void IsingSimImpl::calc_cur_mag()
{
    m_cur_mag = std::accumulate(m_cur_spins.begin(), m_cur_spins.end(), 0);
}
