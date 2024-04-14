#ifndef MISA_KMC_VAC_HASH_H
#define MISA_KMC_VAC_HASH_H

#include <array>
#include "lattice/lattice.h"
#include <functional>
typedef std::function<_type_rate(Lattice *lat_nei, const LatticeTypes::lat_type ghost_atom,
                                 const _type_dir_id _1nn_offset)>
    rateCallback2;

class VacancyHash {
public:

    std::array<double, 8> nn_rates;
    unsigned char avail_trans_dir;

    inline unsigned char availTranDirs2(unsigned char nei_status, Lattice *_1nn_lats);
    void beforeRatesUpdate2(Lattice list_1nn[8], unsigned char status_1nn);

    void updateRates2(Lattice &lattice, Lattice list_1nn[8], unsigned char status_1nn,
                   rateCallback2 callback);
};

#endif // MISA_KMC_VAC_HASH_H