from libc.stdint cimport int8_t, int32_t, uint32_t
from libcpp.vector cimport vector

cdef extern from "ising_impl.hpp":
    cdef cppclass IsingSimImpl:
        void init[S, P, T](vector[S]& neighbors,
            vector[P]& nei_start,
            double beta,
            double J,
            int seed,
            vector[T]& cur_spins) except +

        void iterate(size_t n_iterations)
        void burn_in(size_t n_iterations)
        void iterate[T](size_t n_measures, T aggregator, size_t measure_steps)

        vector[int8_t] get_spins() const
        int32_t get_mag() const
        double get_ene() const
        uint32_t get_N() const
        double get_beta() const
        double get_J() const

