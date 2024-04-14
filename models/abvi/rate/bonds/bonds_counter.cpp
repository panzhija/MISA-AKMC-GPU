//
// Created by genshen on 2019-01-20.
//

#include "bonds_counter.h"
#include <logs/logs.h>
#include <iostream>
/**
 * table 7 in paper {Vincent, E., C. S. Becquart, and C. Domain. \
 * "Solute interaction with point defects in α Fe during thermal ageing: A combined ab initio and atomic kinetic Monte
 * Carlo approach."  \ Journal of nuclear materials 351.1-3 (2006): 88-99}.
 */
const std::map<bonds::PairBond::bond_type, bonds::_type_pair_ia> bonds::BondsCounter::_1nn_bonds = {
    {PairBond::FeFe, -0.778}, {PairBond::VV, 0.315},    {PairBond::CuCu, -0.581}, {PairBond::NiNi, -0.793},
    {PairBond::MnMn, -0.438},

    {PairBond::VFe, -0.161},  {PairBond::FeCu, -0.609}, {PairBond::FeNi, -0.821}, {PairBond::FeMn, -0.648},

    {PairBond::VCu, -0.103},  {PairBond::VNi, -0.234},  {PairBond::VMn, -0.151},

    {PairBond::CuNi, -0.692}, {PairBond::CuMn, -0.519}, {PairBond::NiMn, -0.831},
};

const std::map<bonds::PairBond::bond_type, bonds::_type_pair_ia> bonds::BondsCounter::_2nn_bonds = {
    {PairBond::FeFe, -0.389}, {PairBond::VV, -0.214},   {PairBond::CuCu, -0.389}, {PairBond::NiNi, -0.389},
    {PairBond::MnMn, -0.389},

    {PairBond::VFe, -0.161},  {PairBond::FeCu, -0.344}, {PairBond::FeNi, -0.399}, {PairBond::FeMn, -0.364},

    {PairBond::VCu, -0.206},  {PairBond::VNi, -0.351},  {PairBond::VMn, -0.206},

    {PairBond::CuNi, -0.344}, {PairBond::CuMn, -0.249}, {PairBond::NiMn, -0.464},
};

bonds::_type_pair_ia bonds::BondsCounter::count(LatticesList *lat_list, _type_lattice_id source_id,
                                                LatticeTypes src_atom_type) {
  Lattice _1nn_neighbour[LatticesList::MAX_1NN]; // todo new array many times.
  _type_neighbour_status _1nn_status = lat_list->get1nnStatus(source_id);
  _type_lattice_size x = source_id % lat_list->meta.size_x;
  _type_lattice_size y = (source_id / lat_list->meta.size_x) % lat_list->meta.size_y;
  _type_lattice_size z = source_id / (lat_list->meta.size_x * lat_list->meta.size_y);
  //lat_list->get1nn2(source_id, _1nn_neighbour);
  lat_list->get1nn2(x, y, z, _1nn_neighbour);
  // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
  _type_pair_ia energy = 0.0;
  for (int b = 0; b < LatticesList::MAX_NEI_BITS; b++) {
    if ((_1nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_1nn_neighbour[b].type.isAtom() || _1nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b].type, src_atom_type);
        energy += _1nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 1 is : {} .\n", energy);
  // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
  Lattice _2nn_neighbour[LatticesList::MAX_2NN]; // todo new array many times.
  _type_neighbour_status _2nn_status = lat_list->get2nnStatus(source_id);
  // x = source_id % lat_list->meta.size_x;
  // y = (source_id / lat_list->meta.size_x) % box->lattice_list->meta.size_y;
  // z = source_id / (lat_list->meta.size_x * box->lattice_list->meta.size_y);
  lat_list->get2nn2(x, y, z, _2nn_neighbour);
  // kiwi::logs::v(" ", " x is : {} y is : {} z is : {} .\n", x, y, z);
  // for(int a = 0; a < 6; a++){
  //   kiwi::logs::v(" ", " {} is : {} .\n", a, _2nn_neighbour[a].type._type);
  // }
  // lat_list->get2nn(source_id, _2nn_neighbour);
  for (int b = 0; b < LatticesList::MAX_2NN; b++) {
    if ((_2nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_2nn_neighbour[b].type.isAtom() || _2nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b].type, src_atom_type);
        energy += _2nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 2 is : {} .\n", energy);
  return energy;
}

bonds::_type_pair_ia bonds::BondsCounter::count_exchange(LatticesList *lat_list, _type_lattice_id now_id, _type_lattice_id source_id,
                                                         LatticeTypes atom, LatticeTypes src_atom_type) {
  Lattice _1nn_neighbour[LatticesList::MAX_1NN]; // todo new array many times.
  _type_neighbour_status _1nn_status = lat_list->get1nnStatus(now_id);
  _type_lattice_size x = now_id % lat_list->meta.size_x;
  _type_lattice_size y = (now_id / lat_list->meta.size_x) % lat_list->meta.size_y;
  _type_lattice_size z = now_id / (lat_list->meta.size_x * lat_list->meta.size_y);
  //lat_list->get1nn2(source_id, _1nn_neighbour);
  lat_list->get1nn3(x, y, z, _1nn_neighbour, source_id, atom);
  // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
  _type_pair_ia energy = 0.0;
  for (int b = 0; b < LatticesList::MAX_NEI_BITS; b++) {
    if ((_1nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_1nn_neighbour[b].type.isAtom() || _1nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b].type, src_atom_type);
        energy += _1nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 1 is : {} .\n", energy);
  // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
  Lattice _2nn_neighbour[LatticesList::MAX_2NN]; // todo new array many times.
  _type_neighbour_status _2nn_status = lat_list->get2nnStatus(now_id);
  // x = source_id % lat_list->meta.size_x;
  // y = (source_id / lat_list->meta.size_x) % box->lattice_list->meta.size_y;
  // z = source_id / (lat_list->meta.size_x * box->lattice_list->meta.size_y);
  lat_list->get2nn3(x, y, z, _2nn_neighbour, source_id, atom);
  // kiwi::logs::v(" ", " x is : {} y is : {} z is : {} .\n", x, y, z);
  // for(int a = 0; a < 6; a++){
  //   kiwi::logs::v(" ", " {} is : {} .\n", a, _2nn_neighbour[a].type._type);
  // }
  // lat_list->get2nn(source_id, _2nn_neighbour);
  for (int b = 0; b < LatticesList::MAX_2NN; b++) {
    if ((_2nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_2nn_neighbour[b].type.isAtom() || _2nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b].type, src_atom_type);
        energy += _2nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 2 is : {} .\n", energy);
  return energy;
}