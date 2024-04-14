#include <array>
#include <cassert>
#include "vac_hash.h"
#include <logs/logs.h>
#include <iostream>

unsigned char VacancyHash::availTranDirs2(unsigned char nei_status, Lattice *_1nn_lats){
  unsigned char atom_nei_status = 0;
  for (int b = 0; b < 8; b++) {
    if (((nei_status >> b) & 1) && _1nn_lats[b].type.isAtom()) { // can trans
      atom_nei_status |= 1 << b;
    }
  }
  return atom_nei_status;
}


void VacancyHash::beforeRatesUpdate2(Lattice list_1nn[8], unsigned char status_1nn){
    for (_type_rate &rate : nn_rates) {
      rate = 0;
    }
    // set available transition dir.
    avail_trans_dir = availTranDirs2(status_1nn, list_1nn);
    //kiwi::logs::v(" ", " avail_trans_dir 11111111111 is : {}.\n", avail_trans_dir);
}

void VacancyHash::updateRates2(Lattice &lattice, Lattice list_1nn[8], unsigned char status_1nn,
                          rateCallback2 callback) {
  std::array<double, 8> nn_rates_backup;
  //unsigned char avail_trans_dir_backup = avail_trans_dir;
  for (unsigned char b = 0; b < 8; b++) {
    //kiwi::logs::v(" ", " avail_trans_dir is : {} b is : {}.\n", avail_trans_dir, b);
    if ((avail_trans_dir >> b) & 1) { // the neighbour lattice is available
      nn_rates[b] = callback(&list_1nn[b], list_1nn[b].type._type, b);
      //avail_trans_dir = avail_trans_dir_backup;
    }
  }
  // nn_rates[0] = nn_rates_backup[0];
  // nn_rates[1] = nn_rates_backup[1];
  // nn_rates[2] = nn_rates_backup[2];
  // nn_rates[3] = nn_rates_backup[3];
  // nn_rates[4] = nn_rates_backup[4];
  // nn_rates[5] = nn_rates_backup[5];
  // nn_rates[6] = nn_rates_backup[6];
  // nn_rates[7] = nn_rates_backup[7];
}