#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <type_traits>
#include <cassert>



template<typename T, typename S>
void sort_neighbors(const std::vector<T>& neighbors_start, std::vector<S> neighbors)
{
    const auto nei_begin = neighbors.begin();
    for (uint32_t i = 0, i_end = neighbors_start.size() - 1; i < i_end; ++i) {
        std::sort(nei_begin + neighbors_start[i], nei_begin + neighbors_start[i + 1]);
    }
}
