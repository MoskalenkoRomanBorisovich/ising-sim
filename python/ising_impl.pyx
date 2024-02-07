from ising_impl cimport IsingSimImpl
from libcpp.vector cimport vector

cdef class IsingSim():
    cdef IsingSimImpl _impl
    def __cinit__(self, 
            vector[int]& neighbors,
            vector[int]& nei_start,
            double beta,
            double J,
            int seed,
            vector[int]& cur_spins):

        self._impl.init[int, int, int](
            neighbors,
            nei_start,
            beta,
            J,
            seed,
            cur_spins)

    def iterate(self, size_t n_iterations = 1):
        self._impl.iterate(n_iterations)
    def burn_in(self, size_t n_iterations):
        self._impl.burn_in(n_iterations)

    def get_spins(self):
        return self._impl.get_spins()
    def get_mag(self):
        return self._impl.get_mag()
    def get_ene(self):
        return self._impl.get_ene()
    def get_N(self):
        return self._impl.get_N()
    def get_beta(self):
        return self._impl.get_beta()
    def get_J(self):
        return self._impl.get_J()