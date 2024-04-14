#ifndef DEVICEVACRATESSOVLER_H
#define DEVICEVACRATESSOVLER_H

#include "../src/lattice/lattice_types.h"

struct dev_Vacancy {
    // LatticeTypes type;
    unsigned long id;
    double rates[8];
    int valid; // 改成 1 一个 byte 来存
    // unsigned char avail_trans_dir = 0;
    // dev_vacancy d_vacancy;
};

struct dev_nnLattice {
    LatticeTypes type; // 改成 1 一个 Byte 来存
    unsigned long id;
};

struct dev_tempLattice {
    unsigned long id;
};

struct dev_event {
    unsigned long from_id;
    unsigned long to_id;
    LatticeTypes to_type;
};

#endif /* DEVICEVACRATESSOVLER_H */