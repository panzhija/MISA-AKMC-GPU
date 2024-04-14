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



#ifndef GPU_SIMULATE_H
#define GPU_SIMULATE_H

#define NN_TOTAL 126

#define RE_SCALE_SIZE 2

void firstSectLocalToGpu(_type_lattice_id *vac_idArray_adv, _type_lattice_count vac_adv);

double first_calculate_GPU(_type_lattice_id *vac_idArray_aft, _type_lattice_count vac_aft, _type_lattice_id *vac_idArray_buf_adv, _type_lattice_count vac_adv, int sect);

double calculate_GPU(const comm::Region<comm::_type_lattice_size> region, int sect);

void selectAndPerformEventGPU(double excepted_rand, int rank, _type_lattice_count step, int sect,
                              std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, const unsigned int sector_id,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z);

void calculateNextAdvRegion(int sect);

void sector_final(int sect);

void initialize_gpu(int process_rank);

void gpu_final();

void gpu_prepare(const double v, const double T, LatticesList *p_list, comm::ColoredDomain *_p_domain);

LatticeTypes::lat_type getType(_type_lattice_id latti_id);

void exchange_lattices(dev_event h_event);

void add_ghost(std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
               _type_lattice_id x, _type_lattice_id y, _type_lattice_id z, Lattice lat_to, const unsigned int sector_id);

void initInfo(_type_lattice_count vac_total, _type_lattice_id *vac_idArray,
              dev_Vacancy *h_vacancy, dev_nnLattice *h_nnneighbour);

void add_surface(_type_lattice_id surface_id,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z);

// extern void transSizeToGPU(long meta_size_z, long meta_size_z, long meta_size_z);

#endif /*GPU_SIMULATE_H*/