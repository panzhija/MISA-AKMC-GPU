#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include "gpuError.h"
#include "DeviceVacRatesSolver.h"

#include <cassert>
#include <logs/logs.h>
#include "../models/abvi/kmc.h"
#include "../models/abvi/box.h"
#include "../models/abvi/rate/vacancy_rates_solver.h"
#include "../models/abvi/defect/vacancy_list.h"
#include "../models/abvi/rate/bonds/pair_bond.hpp"
#include "../models/abvi/rate/rates_types.h"
#include "../models/abvi/event.h"

#include "../src/type_define.h"
#include "../src/lattice/lattice.h"
#include "../src/lattice/lattices_list.h"
#include "../src/lattice/lattice_types.h"
#include "../src/lattice/lattice_list_meta.h"
#include <logs/logs.h>
#include <omp.h>

#ifndef GPU_SIMULATE_H
#define GPU_SIMULATE_H

double first_calculate_GPU(LatticesList *lattice_list, _type_lattice_id *vac_idArray, int vac_total);

double calculate_GPU(LatticesList *lattice_list, const comm::Region<comm::_type_lattice_size> region);

void selectAndPerformEventGPU(LatticesList *lattice_list, double excepted_rand, int rank, unsigned long step, int sect,
                              const comm::ColoredDomain *p_domain, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, const unsigned int sector_id);

void sector_final();

void initialize_gpu(int process_rank);

void gpu_final();

void gpu_prepare(const double v, const double T);

LatticeTypes getType(_type_lattice_id latti_id, LatticesList *lattice_list);

void exchange_lattices(LatticesList *lattice_list, _type_lattice_id from_id, _type_lattice_id to_id);

void print_kernel_time();

void add_ghost(const comm::ColoredDomain *p_domain, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
               _type_lattice_id x, _type_lattice_id y, _type_lattice_id z, Lattice lat_to, const unsigned int sector_id);
// extern void transSizeToGPU(long meta_size_z, long meta_size_z, long meta_size_z);

#endif /*GPU_SIMULATE_H*/