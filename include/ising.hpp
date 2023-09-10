#include "ising_impl.hpp"




class IsingSim {
public:
    // IsingSim();


    inline operator bool() const {
        return bool(pimpl);
    };

private:
    std::unique_ptr<IsingSimImpl> pimpl;
}