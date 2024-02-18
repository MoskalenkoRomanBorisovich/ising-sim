from ising_impl cimport IsingSimImpl
from libcpp.vector cimport vector
from mc_lib.observable cimport RealObservable

cdef class IsingSim():
    cdef IsingSimImpl _impl
    cdef RealObservable mag_a
    cdef RealObservable mag_2
    cdef RealObservable mag_4
    cdef RealObservable ene


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

        self.mag_a = RealObservable()
        self.mag_2 = RealObservable()
        self.mag_4 = RealObservable()
        self.ene = RealObservable()
        

    def iterate(self, size_t n_iterations = 1):
        cdef size_t i
        cdef double cur_val
        for i in range(n_iterations):
            self._impl.iterate()
            cur_val = self._impl.get_mag()
            cur_val /= self._impl.get_N()
            self.mag_a.add_measurement(abs(cur_val))
            self.mag_2.add_measurement(cur_val ** 2)
            self.mag_4.add_measurement(cur_val ** 4)
            cur_val = self._impl.get_ene()
            cur_val /= self._impl.get_N()
            self.ene.add_measurement(cur_val)


    def burn_in(self, size_t n_iterations):
        self._impl.burn_in(n_iterations)


    def clear_buf(self):
        self._impl.clear_buf()

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

    @property
    def mag_a (self):
        return self.mag_a
    @property
    def mag_2 (self):
        return self.mag_2
    @property
    def mag_4 (self):
        return self.mag_4
    @property
    def ene (self):
        return self.ene