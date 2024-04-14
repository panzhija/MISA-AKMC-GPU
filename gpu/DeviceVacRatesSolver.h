#ifndef DEVICEVACRATESSOVLER_H
#define DEVICEVACRATESSOVLER_H

#include "../src/lattice/lattice_types.h"
#include "../src/type_define.h"

struct dev_Vacancy {
    // LatticeTypes type;
    _type_lattice_id id;
    double rates[8];
    int valid;
    // unsigned char avail_trans_dir = 0;
    // dev_vacancy d_vacancy;
};

struct dev_nnLattice {
    LatticeTypes type;
    _type_lattice_id id;
};

struct dev_event {
    _type_lattice_id from_id;
    _type_lattice_id to_id;
    LatticeTypes to_type;
};

#endif /* DEVICEVACRATESSOVLER_H */