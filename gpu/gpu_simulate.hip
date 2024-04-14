#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include "gpu_simulate.h"
#include <map>
#include <iostream>

// Global Variables
namespace gpuSimulate {

    // gpu variable
    double *d_sum_rates;
    // long *d_size;
    double *d_config;
    dev_Vacancy *d_vacancy; // 当前子域的所有空位信息
    dev_nnLattice *d_lat_list; // 当前子域所有空位的1nn信息
    dev_nnLattice *d_2nn_neighbour; // 当前子域所有空位的2nn信息
    dev_nnLattice *d_target_lattice_1neighbour; //当前子域所有空位的1nn的1nn信息
    dev_nnLattice *d_target_lattice_2neighbour; //当前子域所有空位的1nn的2nn信息
    // dev_event *d_event;
    // int *d_selected; // 大小为 2 的数组，其第一个元素代表选中的事件中的空位是 h_vacancy_list 数组下的第 d_selected[0] 个元素，它向着它的 d_selected[1] 方向的 1nn跃迁.
    // unsigned long *energy_temp1;
    dev_event *d_selected;
    int *d_event_id;
    // int *d_index;
    // double *d_excepted_rate;
    // LatticeTypes::lat_type *energy_temp2;
    // LatticeTypes::lat_type *type_temp1;
    // LatticeTypes::lat_type *type_temp2;
    // dev_region *d_region;


    // host variable
    int blockSize;
    int numBlocks;
    // _type_lattice_id h_lat_list_id[8];

    double sum_rates;
    int vac_count;
    int vac_valid_count;
    dev_event h_selected;
    int h_event_id;
    // int *h_selected; // 大小为 2 的数组，对应 gpu 端的 d_selected.
    // _type_lattice_id *h_vacancy_list;
    _type_lattice_id to_x;
    _type_lattice_id to_y;
    _type_lattice_id to_z;
    // int last_is_ghost;
    dev_Vacancy h_vacancy_temp;
    dev_nnLattice *h_lat_list_temp;
    // unsigned long *h_lat_list_temp_real;
    dev_nnLattice *h_2nn_neighbour_temp;
    dev_nnLattice *h_target_lattice_1neighbour_temp;
    dev_nnLattice *h_target_lattice_2neighbour_temp;

    double start_time, end_time;
    double total_time1 = 0.0;
    double total_time2 = 0.0;

}

__device__ double _1nnMapToSwitch(bonds::PairBond::bond_type bondType){
    switch (bondType) {
      case bonds::PairBond::FeFe:
        return -0.778;
      case bonds::PairBond::VV:
        return 0.315;
      case bonds::PairBond::CuCu:
        return -0.581;
      case bonds::PairBond::NiNi:
        return -0.793;
      case bonds::PairBond::MnMn:
        return -0.438;
      case bonds::PairBond::VFe:
        return -0.161;
      case bonds::PairBond::FeCu:
        return -0.609;
      case bonds::PairBond::FeNi:
        return -0.821;
      case bonds::PairBond::FeMn:
        return -0.648;
      case bonds::PairBond::VCu:
        return -0.103;
      case bonds::PairBond::VNi:
        return -0.234;
      case bonds::PairBond::VMn:
        return -0.151;
      case bonds::PairBond::CuNi:
        return -0.692;
      case bonds::PairBond::CuMn:
        return -0.519;
      case bonds::PairBond::NiMn:
        return -0.831;
      default:
        // 处理默认情况
        return 0.0;
    }
}

__device__ double _2nnMapToSwitch(bonds::PairBond::bond_type bondType){
    switch (bondType) {
      case bonds::PairBond::FeFe:
        return -0.389;
      case bonds::PairBond::VV:
        return -0.214;
      case bonds::PairBond::CuCu:
        return -0.389;
      case bonds::PairBond::NiNi:
        return -0.389;
      case bonds::PairBond::MnMn:
        return -0.389;
      case bonds::PairBond::VFe:
        return -0.161;
      case bonds::PairBond::FeCu:
        return -0.344;
      case bonds::PairBond::FeNi:
        return -0.399;
      case bonds::PairBond::FeMn:
        return -0.364;
      case bonds::PairBond::VCu:
        return -0.206;
      case bonds::PairBond::VNi:
        return -0.351;
      case bonds::PairBond::VMn:
        return -0.206;
      case bonds::PairBond::CuNi:
        return -0.344;
      case bonds::PairBond::CuMn:
        return -0.249;
      case bonds::PairBond::NiMn:
        return -0.464;
      default:
        // 处理默认情况
        return 0.0;
    }
}

__global__ void updatennLatti(dev_event *d_selected, dev_nnLattice *d_lat_list, 
                              dev_nnLattice *d_2nn_neighbour, dev_nnLattice *d_target_lattice_1neighbour,
                              dev_nnLattice *d_target_lattice_2neighbour,
                              _type_lattice_count vac_total) {
  _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid < vac_total) {
    for(_type_lattice_count i = 0; i < 8; i++){
      if(d_lat_list[tid * 8 + i].id == d_selected->to_id) d_lat_list[tid * 8 + i].type._type = LatticeTypes::V;
      if(d_lat_list[tid * 8 + i].id == d_selected->from_id) d_lat_list[tid * 8 + i].type._type = d_selected->to_type._type;
    }
    for(_type_lattice_count i = 0; i < 6; i++){
      if(d_2nn_neighbour[tid * 6 + i].id == d_selected->to_id) d_2nn_neighbour[tid * 6 + i].type._type = LatticeTypes::V;
      if(d_2nn_neighbour[tid * 6 + i].id == d_selected->from_id) d_2nn_neighbour[tid * 6 + i].type._type = d_selected->to_type._type;
    }
    for(_type_lattice_count i = 0; i < 64; i++){
      if(d_target_lattice_1neighbour[tid * 64 + i].id == d_selected->to_id) d_target_lattice_1neighbour[tid * 64 + i].type._type = LatticeTypes::V;
      if(d_target_lattice_1neighbour[tid * 64 + i].id == d_selected->from_id) d_target_lattice_1neighbour[tid * 64 + i].type._type = d_selected->to_type._type;
    }
    for(_type_lattice_count i = 0; i < 48; i++){
      if(d_target_lattice_2neighbour[tid * 48 + i].id == d_selected->to_id) d_target_lattice_2neighbour[tid * 48 + i].type._type = LatticeTypes::V;
      if(d_target_lattice_2neighbour[tid * 48 + i].id == d_selected->from_id) d_target_lattice_2neighbour[tid * 48 + i].type._type = d_selected->to_type._type;
    }
  }
}

__global__ void calcRatesGPU(dev_Vacancy *d_vacancy, dev_nnLattice *d_lat_list, 
                             dev_nnLattice *d_2nn_neighbour, dev_nnLattice *d_target_lattice_1neighbour,
                             dev_nnLattice *d_target_lattice_2neighbour,
                             double *d_config, double* d_sum_rates, int vac_total) {
    // 获取线程的全局ID
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid == 0){
      *d_sum_rates = 0.0;
    }
    // 判断当前tid是否可用
    // 利用共享内存优化原子加
    __shared__ double tmp_sh;
    __shared__ double attempt_freq;
    __shared__ double temperature;
    //__shared__ double uint16_t big_endian_shift_gpu;
    if (threadIdx.x == 0){
      tmp_sh = 0.0;
      attempt_freq = d_config[0];
      temperature = d_config[1];
      //big_endian_shift_gpu = 8;
    }
    __syncthreads();
    if(tid < vac_total && d_vacancy[tid].valid != 0){
    //if(tid < vac_total){
      // 从tid获取对应的x, y, z索引。
      //unsigned char nei_status = 0xFF;

      d_vacancy[tid].rates[0] = 0;
      d_vacancy[tid].rates[1] = 0;
      d_vacancy[tid].rates[2] = 0;
      d_vacancy[tid].rates[3] = 0;
      d_vacancy[tid].rates[4] = 0;
      d_vacancy[tid].rates[5] = 0;
      d_vacancy[tid].rates[6] = 0;
      d_vacancy[tid].rates[7] = 0;
      double vac_energy = 0;
      int first_ir = 1;
      unsigned char atom_nei_status = 0;
      for (int b = 0; b < 8; b++) {
        if (d_lat_list[tid * 8 + b].type._type <= LatticeTypes::Mn && d_lat_list[tid * 8 + b].type._type != LatticeTypes::V) { // can trans
          atom_nei_status |= 1 << b;
        }
      }

      for (unsigned char b = 0; b < 8; b++) {
        if ((atom_nei_status >> b) & 1) { // the neighbour lattice is available
          double active_energy;
          switch (d_lat_list[tid * 8 + b].type._type) {
            case LatticeTypes::Fe:
              active_energy = 0.62;
              break;
            case LatticeTypes::Cu:
              active_energy = 0.54;
              break;
            case LatticeTypes::Ni:
              active_energy = 0.68;
              break;
            case LatticeTypes::Mn:
              active_energy = 1.03;
              break;
            default:
              active_energy = 0;
              break;
          }
          // calculate system energy before transition.
          // const LatticeTypes::lat_type ghost_atom = d_lat_list[tid * 8 + b].type._type;
          double e_before = 0;
          {
            // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
            double e_src, e_des;
            unsigned char _1nn_status = 0xFF;
            unsigned char _2nn_status = 0x2F;
            double energy = 0.0;
            if(first_ir == 1){
              // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
              for (int b2 = 0; b2 < 8; b2++) {
                if ((_1nn_status >> b2) & 0x1) {
                  // we assume that src_atom_type is single atom or vacancy.
                  if (d_lat_list[tid * 8 + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_1nn = d_lat_list[tid * 8 + b2].type._type;
                    LatticeTypes::lat_type temp_src = LatticeTypes::V;
                    if (temp_1nn > temp_src) {
                      LatticeTypes::lat_type temp = temp_1nn;
                      temp_1nn = temp_src;
                      temp_src = temp;
                    }
                    // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                    const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                    vac_energy += _1nnMapToSwitch(bond);
                  }
                }
              }
          
              // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
              for (int b2 = 0; b2 < 6; b2++) {
                if ((_2nn_status >> b2) & 0x1) {
                  // we assume that src_atom_type is single atom or vacancy.

                  if (d_2nn_neighbour[tid * 6 + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_2nn = d_2nn_neighbour[tid * 6 + b2].type._type;
                    LatticeTypes::lat_type temp_src = LatticeTypes::V;
                    if (temp_2nn > temp_src) {
                      LatticeTypes::lat_type temp = temp_2nn;
                      temp_2nn = temp_src;
                      temp_src = temp;
                    }
                    // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                    const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                    vac_energy += _2nnMapToSwitch(bond);
                  }
                }
              }
              first_ir = 0;
            }
            e_src = vac_energy;
//_______________________________________
            energy = 0.0;
            // double e_src = bonds::BondsCounter::count(&lattice_list, source_id, source_lattice.type);
            // bonds energy of des lattice contributed by its 1nn/2nn neighbour lattice(it is an atom).
              
            // Lattice *_1nn_neighbour2[LatticesList::MAX_1NN] = {nullptr}; // todo new array many times.
            // long x2, y2, z2;
            //nunsigned long lid = target_lat_id;
            // x3 = lid2 % meta_size_x;
            // lid2 = lid2 / meta_size_x;
            // y3 = lid2 % meta_size_y;
            // z3 = lid2 / meta_size_y;
                // meta.getCoordByLId(source_lat_id, &x3, &y3, &z3);
                // return get1nnStatus(x, y, z);
              
            // _type_neighbour_status _1nn_status = 0xFF;
            // lat_list->get1nn(source_id, _1nn_neighbour);

            // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_lat_list[tid * 8 + b].type._type;
                  if (temp_1nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_1nn;
                    temp_1nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                  energy += _1nnMapToSwitch(bond);
                }
              }
            }
            // if(is_first) *energy_temp1 = energy;
            // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
            // Lattice *_2nn_neighbour[LatticesList::MAX_2NN] = {nullptr}; // todo new array many times.
            // _type_neighbour_status _2nn_status = 0x2F;
                // lat_list->get2nn(source_id, _2nn_neighbour);
            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_target_lattice_2neighbour[tid * 48 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_target_lattice_2neighbour[tid * 48 + b * 6 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_lat_list[tid * 8 + b].type._type;
                  if (temp_2nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_2nn;
                    temp_2nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                  energy += _2nnMapToSwitch(bond);
                }
              }
            }
            // if(is_first) *energy_temp2 = energy;
            // is_first = false;
            e_des = energy;

            // double e_des = bonds::BondsCounter::count(&lattice_list, target_id, target_lattice.type);
            e_before = e_src + e_des;
            // energy_temp2[b] = energy;
          }
          // exchange atoms of vacancy and neighbour lattice atom/vacancy.
          // note: ghost_atom equals to origin target_lattice.type
          // lattice->type._type = ghost_atom; // const LatticeTypes::lat_type ghost_atom = lat_list[b]->type._type;
          // lat_list[b]->type._type = LatticeTypes::V;

// ______________________________________________________________________________________________________________________________________



          // calculate system energy after transition.

          double e_after = 0;
          {
            // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
            double e_src, e_des;
        
            // _type_lattice_coord x3, y3, z3;
            // _type_lattice_id lid2;
            // lid2 = source_lat_id;
            // x3 = lid2 % d_vac_rate->lattice_list.meta.size_x;
            // lid2 = lid2 / d_vac_rate->lattice_list.meta.size_x;
            // y3 = lid2 % d_vac_rate->lattice_list.meta.size_y;
            // z3 = lid2 / d_vac_rate->lattice_list.meta.size_y;
            // meta.getCoordByLId(source_lat_id, &x3, &y3, &z3);
            // return get1nnStatus(x, y, z);    
            unsigned char _1nn_status = 0xFF;
            // lat_list->get1nn(source_id, _1nn_neighbour);

            // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
                
            double energy = 0.0;

            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                LatticeTypes::lat_type temp_type = d_lat_list[tid * 8 + b2].type._type;
                if (d_lat_list[tid * 8 + b2].id == d_lat_list[tid * 8 + b].id){
                  temp_type = LatticeTypes::V;
                }
                // if (b2 == b){
                //   temp_type = LatticeTypes::V;
                // }
                if (temp_type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = temp_type; 
                  LatticeTypes::lat_type temp_src = d_lat_list[tid * 8 + b].type._type;
                  if (temp_1nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_1nn;
                    temp_1nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                  bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                  energy += _1nnMapToSwitch(bond);
                }
              }
            }

            // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
            unsigned char _2nn_status = 0x2F;
            
            // lat_list->get2nn(source_id, _2nn_neighbour);
            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_2nn_neighbour[tid * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_2nn_neighbour[tid * 6 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_lat_list[tid * 8 + b].type._type;
                  if (temp_2nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_2nn;
                    temp_2nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                  energy += _2nnMapToSwitch(bond);
                }
              }
            }
            e_src = energy;
//_______________________________________
            energy = 0.0;
            // double e_src = bonds::BondsCounter::count(&lattice_list, source_id, source_lattice.type);
            // bonds energy of des lattice contributed by its 1nn/2nn neighbour lattice(it is an atom).
              
            // Lattice *_1nn_neighbour2[LatticesList::MAX_1NN] = {nullptr}; // todo new array many times.
            // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                LatticeTypes::lat_type temp_type = d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].type._type;
                if (d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].id == d_vacancy[tid].id){
                  temp_type = d_lat_list[tid * 8 + b].type._type;
                }
                // if (b2 + b == 7) {
                //   temp_type = d_lat_list[tid * 8 + b].type._type;
                // }
                // if (d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].id == d_lattices[tid].id){
                //   temp_type = d_lat_list[tid * 8 + b].type._type;
                // }
                if (temp_type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = temp_type;
                  LatticeTypes::lat_type temp_src = LatticeTypes::V;
                  if (temp_1nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_1nn;
                    temp_1nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                  energy += _1nnMapToSwitch(bond);
                }
              }
            }
            // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
            // Lattice *_2nn_neighbour[LatticesList::MAX_2NN] = {nullptr}; // todo new array many times.
            // _type_neighbour_status _2nn_status = 0x2F;
            // lat_list->get2nn(source_id, _2nn_neighbour);
            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_target_lattice_2neighbour[tid * 48 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_target_lattice_2neighbour[tid * 48 + b * 6 + b2].type._type;;
                  LatticeTypes::lat_type temp_src = LatticeTypes::V;
                  if (temp_2nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_2nn;
                    temp_2nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                  energy += _2nnMapToSwitch(bond);
                }
              }
            }
            e_des = energy;
            // double e_des = bonds::BondsCounter::count(&lattice_list, target_id, target_lattice.type);
            e_after = e_src + e_des;
          }
          active_energy += (e_after - e_before) / 2;
          double rate = attempt_freq * expf(-active_energy / (8.617330350E-5 * temperature));
          d_vacancy[tid].rates[b] = rate;
          // return arrhenius(attempt_freq, temperature, active_energy);

          // _type_rate rate = d_vac_rate.rate(lattice, *lat_nei, ghost_atom, _1nn_offset);
          atomicAdd(&tmp_sh, rate);
        }
      }
    }
  __syncthreads();
  if (threadIdx.x == 0) atomicAdd(d_sum_rates, tmp_sh);
}


__global__ void selectEventGPU(dev_Vacancy *d_vacancy, dev_nnLattice *d_lat_list, dev_event *d_selected, 
                               int *d_event_id, int vac_count, double *d_sum_rates, double excepted_rand) {

  // 获取线程的全局ID
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid == 0){
    __shared__ double excepted_rate;
    __shared__ double rate_accumulator;
    excepted_rate = *d_sum_rates * excepted_rand;
    rate_accumulator = 0;

    // 线性搜索合适的事件 平均时间复杂度为 O( N/2 + 4) N为空位数量
    int tag = 1;
    for(int i = 0; i < vac_count; i++){
      if(d_vacancy[i].valid != 0) {
        if(tag == 1){
          for (int rate_index = 0; rate_index < 8; rate_index++) {
            rate_accumulator += d_vacancy[i].rates[rate_index];
            if(rate_accumulator > excepted_rate) {
              if(tag == 1){
                d_selected->from_id = d_vacancy[i].id;
                d_selected->to_id = d_lat_list[i * 8 + rate_index].id;
                d_selected->to_type._type = d_lat_list[i * 8 + rate_index].type._type;
                *d_event_id = i;
                // *d_index = rate_index;
                // *d_excepted_rate = excepted_rate;
                tag = 0;
                break;
              }
            }
          }
        }else{
          break;
        }
      }
    }
  }
}

double first_calculate_GPU(LatticesList *lattice_list, _type_lattice_id *vac_idArray, int vac_total){
  
  using namespace gpuSimulate;
    
    // size_t free_mem, total_mem;
    // HANDLE_HIP(hipMemGetInfo(&free_mem, &total_mem));
    // float free_memory = free_mem / (1024.0 * 1024.0);
    // float total_memory = total_mem / (1024.0 * 1024.0);
    // kiwi::logs::v(" ", "Device_mem_before: free = {} total = {}.\n", free_memory, total_memory);
    // for(int i = 0; i < vac_total; i++){
    //   kiwi::logs::v(" ", " id is : {} .\n", vac_idArray[i]);
    // }
    sum_rates = 0;
    vac_count = vac_total;
    vac_valid_count = vac_total;
    
    // h_vacancy_list = (_type_lattice_id *)malloc(vac_total * sizeof(_type_lattice_id));

    dev_Vacancy *h_vacancy;
    dev_nnLattice *h_lat_list;
    // unsigned long *h_lat_list_real;
    dev_nnLattice *h_2nn_neighbour;
    dev_nnLattice *h_target_lattice_1neighbour;
    dev_nnLattice *h_target_lattice_2neighbour;

    h_vacancy = (dev_Vacancy *)malloc(vac_total * sizeof(dev_Vacancy));
    h_lat_list = (dev_nnLattice *)malloc(8 * vac_total * sizeof(dev_nnLattice));
    // h_lat_list_real = (unsigned long *)malloc(8 * vac_total * sizeof(unsigned long));
    h_2nn_neighbour = (dev_nnLattice *)malloc(6 * vac_total * sizeof(dev_nnLattice));
    h_target_lattice_1neighbour = (dev_nnLattice *)malloc(64 * vac_total * sizeof(dev_nnLattice));
    h_target_lattice_2neighbour = (dev_nnLattice *)malloc(48 * vac_total * sizeof(dev_nnLattice));

    HANDLE_HIP(hipMalloc((void**)&d_vacancy, vac_total * sizeof(dev_Vacancy)));
    HANDLE_HIP(hipMalloc((void**)&d_lat_list, 8 * vac_total * sizeof(dev_nnLattice)));
    // HANDLE_HIP(hipMalloc((void**)&d_lat_list_real, 8 * vac_total * sizeof(unsigned long)));
    HANDLE_HIP(hipMalloc((void**)&d_2nn_neighbour, 6 * vac_total * sizeof(dev_nnLattice)));
    HANDLE_HIP(hipMalloc((void**)&d_target_lattice_1neighbour, 64 * vac_total * sizeof(dev_nnLattice)));
    HANDLE_HIP(hipMalloc((void**)&d_target_lattice_2neighbour, 48 * vac_total * sizeof(dev_nnLattice)));

    // _type_lattice_id latti_id;

    // 准备可能用到的lattice
    //#pragma omp parallel for num_threads(8)
    for(int i = 0; i < vac_total; i++){
      //if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", "  omp_get_thread_num() is : {} . \n", omp_get_thread_num());
      _type_lattice_id x = vac_idArray[i] % lattice_list->meta.size_x;
      _type_lattice_id y = (vac_idArray[i] / lattice_list->meta.size_x) % lattice_list->meta.size_y;
      _type_lattice_id z = vac_idArray[i] / (lattice_list->meta.size_x * lattice_list->meta.size_y);

      h_vacancy[i].id = vac_idArray[i];
      h_vacancy[i].valid = 1;

      if(x % 2 == 0){

        h_lat_list[i * 8].id = lattice_list->getId(x - 1, y - 1, z - 1);
        h_lat_list[i * 8].type = getType(h_lat_list[i * 8].id, lattice_list);

        h_lat_list[i * 8 + 1].id = lattice_list->getId(x - 1, y - 1, z);
        h_lat_list[i * 8 + 1].type = getType(h_lat_list[i * 8 + 1].id, lattice_list);

        h_lat_list[i * 8 + 2].id = lattice_list->getId(x - 1, y, z - 1);
        h_lat_list[i * 8 + 2].type = getType(h_lat_list[i * 8 + 2].id, lattice_list);

        h_lat_list[i * 8 + 3].id = lattice_list->getId(x - 1, y, z);
        h_lat_list[i * 8 + 3].type = getType(h_lat_list[i * 8 + 3].id, lattice_list);

        h_lat_list[i * 8 + 4].id = lattice_list->getId(x + 1, y - 1, z - 1);
        h_lat_list[i * 8 + 4].type = getType(h_lat_list[i * 8 + 4].id, lattice_list);

        h_lat_list[i * 8 + 5].id = lattice_list->getId(x + 1, y - 1, z);
        h_lat_list[i * 8 + 5].type = getType(h_lat_list[i * 8 + 5].id, lattice_list);

        h_lat_list[i * 8 + 6].id = lattice_list->getId(x + 1, y, z - 1);
        h_lat_list[i * 8 + 6].type = getType(h_lat_list[i * 8 + 6].id, lattice_list);

        h_lat_list[i * 8 + 7].id = lattice_list->getId(x + 1, y, z);
        h_lat_list[i * 8 + 7].type = getType(h_lat_list[i * 8 + 7].id, lattice_list);

      }else{
        
        h_lat_list[i * 8].id = lattice_list->getId(x - 1, y, z);
        h_lat_list[i * 8].type = getType(h_lat_list[i * 8].id, lattice_list);

        h_lat_list[i * 8 + 1].id = lattice_list->getId(x - 1, y, z + 1);
        h_lat_list[i * 8 + 1].type = getType(h_lat_list[i * 8 + 1].id, lattice_list);

        h_lat_list[i * 8 + 2].id = lattice_list->getId(x - 1, y + 1, z);
        h_lat_list[i * 8 + 2].type = getType(h_lat_list[i * 8 + 2].id, lattice_list);

        h_lat_list[i * 8 + 3].id = lattice_list->getId(x - 1, y + 1, z + 1);
        h_lat_list[i * 8 + 3].type = getType(h_lat_list[i * 8 + 3].id, lattice_list);

        h_lat_list[i * 8 + 4].id = lattice_list->getId(x + 1, y, z);
        h_lat_list[i * 8 + 4].type = getType(h_lat_list[i * 8 + 4].id, lattice_list);

        h_lat_list[i * 8 + 5].id = lattice_list->getId(x + 1, y, z + 1);
        h_lat_list[i * 8 + 5].type = getType(h_lat_list[i * 8 + 5].id, lattice_list);

        h_lat_list[i * 8 + 6].id = lattice_list->getId(x + 1, y + 1, z);
        h_lat_list[i * 8 + 6].type = getType(h_lat_list[i * 8 + 6].id, lattice_list);

        h_lat_list[i * 8 + 7].id = lattice_list->getId(x + 1, y + 1, z + 1);
        h_lat_list[i * 8 + 7].type = getType(h_lat_list[i * 8 + 7].id, lattice_list);

      }

      h_2nn_neighbour[i * 6].id = lattice_list->getId(x - 2, y, z);
      h_2nn_neighbour[i * 6].type = getType(h_2nn_neighbour[i * 6].id, lattice_list);

      h_2nn_neighbour[i * 6 + 1].id = lattice_list->getId(x, y - 1, z);
      h_2nn_neighbour[i * 6 + 1].type = getType(h_2nn_neighbour[i * 6 + 1].id, lattice_list);

      h_2nn_neighbour[i * 6 + 2].id = lattice_list->getId(x, y, z - 1);
      h_2nn_neighbour[i * 6 + 2].type = getType(h_2nn_neighbour[i * 6 + 2].id, lattice_list);

      h_2nn_neighbour[i * 6 + 3].id = lattice_list->getId(x, y, z + 1);
      h_2nn_neighbour[i * 6 + 3].type = getType(h_2nn_neighbour[i * 6 + 3].id, lattice_list);

      h_2nn_neighbour[i * 6 + 4].id = lattice_list->getId(x, y + 1, z);
      h_2nn_neighbour[i * 6 + 4].type = getType(h_2nn_neighbour[i * 6 + 4].id, lattice_list);

      h_2nn_neighbour[i * 6 + 5].id = lattice_list->getId(x + 2, y, z);
      h_2nn_neighbour[i * 6 + 5].type = getType(h_2nn_neighbour[i * 6 + 5].id, lattice_list);

      for(int b = 0; b < 8; b++){
        _type_lattice_id x2 = h_lat_list[i * 8 + b].id % lattice_list->meta.size_x;
        _type_lattice_id y2 = (h_lat_list[i * 8 + b].id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
        _type_lattice_id z2 = h_lat_list[i * 8 + b].id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

        if(x2 % 2 == 0){
          
          h_target_lattice_1neighbour[i * 64 + b * 8].id = lattice_list->getId(x2 - 1, y2 - 1, z2 - 1);
          h_target_lattice_1neighbour[i * 64 + b * 8].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 1].id = lattice_list->getId(x2 - 1, y2 - 1, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 1].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 1].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 2].id = lattice_list->getId(x2 - 1, y2, z2 - 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 2].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 2].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 3].id = lattice_list->getId(x2 - 1, y2, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 3].type = getType( h_target_lattice_1neighbour[i * 64 + b * 8 + 3].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 4].id = lattice_list->getId(x2 + 1, y2 - 1, z2 - 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 4].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 4].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 5].id = lattice_list->getId(x2 + 1, y2 - 1, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 5].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 5].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 6].id = lattice_list->getId(x2 + 1, y2, z2 - 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 6].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 6].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 7].id = lattice_list->getId(x2 + 1, y2, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 7].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 7].id, lattice_list);

        }else{

          h_target_lattice_1neighbour[i * 64 + b * 8].id = lattice_list->getId(x2 - 1, y2, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 1].id = lattice_list->getId(x2 - 1, y2, z2 + 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 1].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 1].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 2].id = lattice_list->getId(x2 - 1, y2 + 1, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 2].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 2].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 3].id = lattice_list->getId(x2 - 1, y2 + 1, z2 + 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 3].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 3].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 4].id = lattice_list->getId(x2 + 1, y2, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 4].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 4].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 5].id = lattice_list->getId(x2 + 1, y2, z2 + 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 5].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 5].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 6].id = lattice_list->getId(x2 + 1, y2 + 1, z2);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 6].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 6].id, lattice_list);

          h_target_lattice_1neighbour[i * 64 + b * 8 + 7].id = lattice_list->getId(x2 + 1, y2 + 1, z2 + 1);
          h_target_lattice_1neighbour[i * 64 + b * 8 + 7].type = getType(h_target_lattice_1neighbour[i * 64 + b * 8 + 7].id, lattice_list);

        }

        h_target_lattice_2neighbour[i * 48 + b * 6].id = lattice_list->getId(x2 - 2, y2, z2);
        h_target_lattice_2neighbour[i * 48 + b * 6].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6].id, lattice_list);

        h_target_lattice_2neighbour[i * 48 + b * 6 + 1].id = lattice_list->getId(x2, y2 - 1, z2);
        h_target_lattice_2neighbour[i * 48 + b * 6 + 1].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6 + 1].id, lattice_list);

        h_target_lattice_2neighbour[i * 48 + b * 6 + 2].id = lattice_list->getId(x2, y2, z2 - 1);
        h_target_lattice_2neighbour[i * 48 + b * 6 + 2].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6 + 2].id, lattice_list);

        h_target_lattice_2neighbour[i * 48 + b * 6 + 3].id = lattice_list->getId(x2, y2, z2 + 1);
        h_target_lattice_2neighbour[i * 48 + b * 6 + 3].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6 + 3].id, lattice_list);

        h_target_lattice_2neighbour[i * 48 + b * 6 + 4].id = lattice_list->getId(x2, y2 + 1, z2);
        h_target_lattice_2neighbour[i * 48 + b * 6 + 4].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6 + 4].id, lattice_list);

        h_target_lattice_2neighbour[i * 48 + b * 6 + 5].id = lattice_list->getId(x2 + 2, y2, z2);
        h_target_lattice_2neighbour[i * 48 + b * 6 + 5].type = getType(h_target_lattice_2neighbour[i * 48 + b * 6 + 5].id, lattice_list);

      }
    }

    HANDLE_HIP(hipMemcpy(d_vacancy, h_vacancy, vac_total * sizeof(dev_Vacancy), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_lat_list, h_lat_list, 8 * vac_total * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    // HANDLE_HIP(hipMemcpy(d_lat_list_real, h_lat_list_real, 8 * vac_total * sizeof(unsigned long), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_2nn_neighbour, h_2nn_neighbour, 6 * vac_total * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_target_lattice_1neighbour, h_target_lattice_1neighbour, 64 * vac_total * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_target_lattice_2neighbour, h_target_lattice_2neighbour, 48 * vac_total * sizeof(dev_nnLattice), hipMemcpyHostToDevice));

    blockSize = 64;
    numBlocks = (vac_count + blockSize - 1) / blockSize;


    //start_time = MPI_Wtime();
    // 调用GPU内核函数
    hipLaunchKernelGGL(calcRatesGPU, numBlocks, blockSize, 0, 0, 
                  d_vacancy,
                  d_lat_list,
                  d_2nn_neighbour, 
                  d_target_lattice_1neighbour,
                  d_target_lattice_2neighbour, 
                  d_config, 
                  d_sum_rates, 
                  vac_count);
    // hipDeviceSynchronize();
    // end_time = MPI_Wtime() - start_time;
    // total_time1 += end_time;
    CHECK_FOR_ERROR();

    HANDLE_HIP(hipMemcpy(&sum_rates, d_sum_rates, sizeof(double), hipMemcpyDeviceToHost));

    // for (int aaa = 0; aaa < vac_count; aaa++) {
    //   std::cout << "  " << h_lattice[aaa].id << "  ";
    // }
    // std::cout << "" << std::endl;

    // HANDLE_HIP(hipMemGetInfo(&free_mem, &total_mem));
    // free_memory = free_mem / (1024.0 * 1024.0);
    // total_memory = total_mem / (1024.0 * 1024.0);
    // kiwi::logs::v(" ", "Device_mem_after: free = {} total = {}.\n", free_memory, total_memory);

    free(h_vacancy);
    free(h_lat_list);
    // free(h_lat_list_real);
    free(h_2nn_neighbour);
    free(h_target_lattice_1neighbour);
    free(h_target_lattice_2neighbour);
    free(vac_idArray);
    return sum_rates;
}

void print_kernel_time() {
  using namespace gpuSimulate;
  kiwi::logs::v(" ", " kernel time is : {}.\n", total_time1);
}

void selectAndPerformEventGPU(LatticesList *lattice_list, double excepted_rand, int rank, unsigned long step, int sect,
                              const comm::ColoredDomain *p_domain, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, const unsigned int sector_id){
  using namespace gpuSimulate;

  blockSize = 1;
  numBlocks = 1;

  // dev_Vacancy *h_vacancy_temp;
  // h_vacancy_temp = (dev_Vacancy *)malloc(vac_count * sizeof(dev_Vacancy));
  // HANDLE_HIP(hipMemcpy(h_vacancy_temp, d_vacancy, vac_count * sizeof(dev_Vacancy), hipMemcpyDeviceToHost));
  // for(int i = 0; i < vac_count; i++){
  //   if(h_vacancy_temp[i].valid == 1){
  //     for(int b = 0; b < 8; b++){
  //       kiwi::logs::v(" ", " id is : {} rates {} is {} .\n", h_vacancy_temp[i].id, b, h_vacancy_temp[i].rates[b]);
  //     }
  //   }
  // }
  // free(h_vacancy_temp);

  hipLaunchKernelGGL(selectEventGPU, numBlocks, blockSize, 0, 0, 
                  d_vacancy,
                  d_lat_list,
                  d_selected,
                  d_event_id,
                  vac_count,
                  d_sum_rates,
                  excepted_rand);
  
  CHECK_FOR_ERROR();

  //kiwi::logs::v("", "excepted_rate is : {} .\n", sum_rates * excepted_rand);

  // 有事件被选中，则交换两个lattice的type

  HANDLE_HIP(hipMemcpy(&h_selected, d_selected, sizeof(dev_event), hipMemcpyDeviceToHost));

  // int index;
  // double h_excepted_rate;
  // HANDLE_HIP(hipMemcpy(&index, d_index, sizeof(int), hipMemcpyDeviceToHost));
  // HANDLE_HIP(hipMemcpy(&h_excepted_rate, d_excepted_rate, sizeof(double), hipMemcpyDeviceToHost));
  // kiwi::logs::v(" ", " excepted_rate is : {} selected index is : {} .\n", h_excepted_rate, index);
  // last_is_ghost = 0;
  _type_lattice_size from_x = h_selected.from_id % lattice_list->meta.size_x;
  _type_lattice_size from_y = (h_selected.from_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  _type_lattice_size from_z = h_selected.from_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

  to_x = h_selected.to_id % lattice_list->meta.size_x;
  to_y = (h_selected.to_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  to_z = h_selected.to_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);
  exchange_lattices(lattice_list, h_selected.from_id, h_selected.to_id);
  // h_vacancy_list[h_selected[0]] = to_id;
  if(lattice_list->meta.isGhostLat(h_selected.to_id)){
    // last_is_ghost = 1;
    Lattice lat_to;
    lat_to.id = h_selected.to_id;
    lat_to.type = LatticeTypes{LatticeTypes::V};
    add_ghost(p_domain, exchange_ghost, to_x, to_y, to_z, lat_to, sector_id);
  }
  // kiwi::logs::v(" ", " from_id is : {} to_id is : {} .\n", h_selected.from_id, h_selected.to_id);
}

double calculate_GPU(LatticesList *lattice_list, const comm::Region<comm::_type_lattice_size> region){
  using namespace gpuSimulate;

  HANDLE_HIP(hipMemcpy(&h_event_id, d_event_id, sizeof(int), hipMemcpyDeviceToHost));

  _type_lattice_size x3 = to_x - lattice_list->meta.ghost_x;
  _type_lattice_size y3 = to_y - lattice_list->meta.ghost_y;
  _type_lattice_size z3 = to_z - lattice_list->meta.ghost_z;

  // 判断转移后的位置是否位于阴影区 或者 是否还位于当前子域内。
  if (x3 < 0 || y3 < 0 || z3 < 0 || x3 >= lattice_list->meta.box_x || y3 >= lattice_list->meta.box_y || z3 >= lattice_list->meta.box_z || !(2 * region.x_low <= to_x && to_x < 2 * region.x_high && region.y_low <= to_y && to_y < region.y_high && region.z_low <= to_z && to_z < region.z_high)){
    // 转移后的位置位于阴影区 或者 位于其他子域内。
    vac_valid_count--;
    if(vac_valid_count > 0){
      // 当前子域还有空位
      h_vacancy_temp.valid = 0;
      HANDLE_HIP(hipMemcpy(d_vacancy + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
    }else{
      // 当前子域无空位
      return 0;
    }
  
  }else{
    // 转移后的位置没有位于阴影区 并且 还在当前子域内。
    // 更新相关数据。

    h_vacancy_temp.id = h_selected.to_id;
    h_vacancy_temp.valid = 1;

    // _type_lattice_id latti_id;

    if(to_x % 2 == 0){
      
      h_lat_list_temp[0].id = lattice_list->getId(to_x - 1, to_y - 1, to_z - 1);
      h_lat_list_temp[0].type = getType(h_lat_list_temp[0].id, lattice_list);

      h_lat_list_temp[1].id = lattice_list->getId(to_x - 1, to_y - 1, to_z);
      h_lat_list_temp[1].type = getType(h_lat_list_temp[1].id, lattice_list);

      h_lat_list_temp[2].id = lattice_list->getId(to_x - 1, to_y, to_z - 1);
      h_lat_list_temp[2].type = getType(h_lat_list_temp[2].id, lattice_list);

      h_lat_list_temp[3].id = lattice_list->getId(to_x - 1, to_y, to_z);
      h_lat_list_temp[3].type = getType(h_lat_list_temp[3].id, lattice_list);

      h_lat_list_temp[4].id = lattice_list->getId(to_x + 1, to_y - 1, to_z - 1);
      h_lat_list_temp[4].type = getType(h_lat_list_temp[4].id, lattice_list);

      h_lat_list_temp[5].id = lattice_list->getId(to_x + 1, to_y - 1, to_z);
      h_lat_list_temp[5].type = getType(h_lat_list_temp[5].id, lattice_list);

      h_lat_list_temp[6].id = lattice_list->getId(to_x + 1, to_y, to_z - 1);
      h_lat_list_temp[6].type = getType(h_lat_list_temp[6].id, lattice_list);

      h_lat_list_temp[7].id = lattice_list->getId(to_x + 1, to_y, to_z);
      h_lat_list_temp[7].type = getType(h_lat_list_temp[7].id, lattice_list);

    }else{

      h_lat_list_temp[0].id = lattice_list->getId(to_x - 1, to_y, to_z);
      h_lat_list_temp[0].type = getType(h_lat_list_temp[0].id, lattice_list);

      h_lat_list_temp[1].id = lattice_list->getId(to_x - 1, to_y, to_z + 1);
      h_lat_list_temp[1].type = getType(h_lat_list_temp[1].id, lattice_list);

      h_lat_list_temp[2].id = lattice_list->getId(to_x - 1, to_y + 1, to_z);
      h_lat_list_temp[2].type = getType(h_lat_list_temp[2].id, lattice_list);

      h_lat_list_temp[3].id = lattice_list->getId(to_x - 1, to_y + 1, to_z + 1);
      h_lat_list_temp[3].type = getType(h_lat_list_temp[3].id, lattice_list);

      h_lat_list_temp[4].id = lattice_list->getId(to_x + 1, to_y, to_z);
      h_lat_list_temp[4].type = getType(h_lat_list_temp[4].id, lattice_list);

      h_lat_list_temp[5].id = lattice_list->getId(to_x + 1, to_y, to_z + 1);
      h_lat_list_temp[5].type = getType(h_lat_list_temp[5].id, lattice_list);

      h_lat_list_temp[6].id = lattice_list->getId(to_x + 1, to_y + 1, to_z);
      h_lat_list_temp[6].type = getType(h_lat_list_temp[6].id, lattice_list);

      h_lat_list_temp[7].id = lattice_list->getId(to_x + 1, to_y + 1, to_z + 1);
      h_lat_list_temp[7].type = getType(h_lat_list_temp[7].id, lattice_list);

    }

    h_2nn_neighbour_temp[0].id = lattice_list->getId(to_x - 2, to_y, to_z);
    h_2nn_neighbour_temp[0].type = getType(h_2nn_neighbour_temp[0].id, lattice_list);

    h_2nn_neighbour_temp[1].id = lattice_list->getId(to_x, to_y - 1, to_z);
    h_2nn_neighbour_temp[1].type = getType(h_2nn_neighbour_temp[1].id, lattice_list);

    h_2nn_neighbour_temp[2].id = lattice_list->getId(to_x, to_y, to_z - 1);
    h_2nn_neighbour_temp[2].type = getType(h_2nn_neighbour_temp[2].id, lattice_list);

    h_2nn_neighbour_temp[3].id = lattice_list->getId(to_x, to_y, to_z + 1);
    h_2nn_neighbour_temp[3].type = getType(h_2nn_neighbour_temp[3].id, lattice_list);

    h_2nn_neighbour_temp[4].id = lattice_list->getId(to_x, to_y + 1, to_z);
    h_2nn_neighbour_temp[4].type = getType(h_2nn_neighbour_temp[4].id, lattice_list);

    h_2nn_neighbour_temp[5].id = lattice_list->getId(to_x + 2, to_y, to_z);
    h_2nn_neighbour_temp[5].type = getType(h_2nn_neighbour_temp[5].id, lattice_list);
        
    //#pragma omp parallel for num_threads(8)
    //start_time = MPI_Wtime();
    for(int b = 0; b < 8; b++){
      _type_lattice_id x4 = h_lat_list_temp[b].id % lattice_list->meta.size_x;
      _type_lattice_id y4 = (h_lat_list_temp[b].id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
      _type_lattice_id z4 = h_lat_list_temp[b].id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

      if(x4 % 2 == 0){
        
        h_target_lattice_1neighbour_temp[b * 8].id = lattice_list->getId(x4 - 1, y4 - 1, z4 - 1);
        h_target_lattice_1neighbour_temp[b * 8].type = getType(h_target_lattice_1neighbour_temp[b * 8].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 1].id = lattice_list->getId(x4 - 1, y4 - 1, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 1].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 1].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 2].id = lattice_list->getId(x4 - 1, y4, z4 - 1);
        h_target_lattice_1neighbour_temp[b * 8 + 2].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 2].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 3].id = lattice_list->getId(x4 - 1, y4, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 3].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 3].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 4].id = lattice_list->getId(x4 + 1, y4 - 1, z4 - 1);
        h_target_lattice_1neighbour_temp[b * 8 + 4].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 4].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 5].id = lattice_list->getId(x4 + 1, y4 - 1, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 5].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 5].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 6].id = lattice_list->getId(x4 + 1, y4, z4 - 1);
        h_target_lattice_1neighbour_temp[b * 8 + 6].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 6].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 7].id = lattice_list->getId(x4 + 1, y4, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 7].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 7].id, lattice_list);
          
      }else{

        h_target_lattice_1neighbour_temp[b * 8].id = lattice_list->getId(x4 - 1, y4, z4);
        h_target_lattice_1neighbour_temp[b * 8].type = getType(h_target_lattice_1neighbour_temp[b * 8].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 1].id = lattice_list->getId(x4 - 1, y4, z4 + 1);
        h_target_lattice_1neighbour_temp[b * 8 + 1].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 1].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 2].id = lattice_list->getId(x4 - 1, y4 + 1, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 2].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 2].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 3].id = lattice_list->getId(x4 - 1, y4 + 1, z4 + 1);
        h_target_lattice_1neighbour_temp[b * 8 + 3].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 3].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 4].id = lattice_list->getId(x4 + 1, y4, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 4].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 4].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 5].id = lattice_list->getId(x4 + 1, y4, z4 + 1);
        h_target_lattice_1neighbour_temp[b * 8 + 5].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 5].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 6].id = lattice_list->getId(x4 + 1, y4 + 1, z4);
        h_target_lattice_1neighbour_temp[b * 8 + 6].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 6].id, lattice_list);

        h_target_lattice_1neighbour_temp[b * 8 + 7].id = lattice_list->getId(x4 + 1, y4 + 1, z4 + 1);
        h_target_lattice_1neighbour_temp[b * 8 + 7].type = getType(h_target_lattice_1neighbour_temp[b * 8 + 7].id, lattice_list);

      }

        h_target_lattice_2neighbour_temp[b * 6].id = lattice_list->getId(x4 - 2, y4, z4);
        h_target_lattice_2neighbour_temp[b * 6].type = getType(h_target_lattice_2neighbour_temp[b * 6].id, lattice_list);

        h_target_lattice_2neighbour_temp[b * 6 + 1].id = lattice_list->getId(x4, y4 - 1, z4);
        h_target_lattice_2neighbour_temp[b * 6 + 1].type = getType(h_target_lattice_2neighbour_temp[b * 6 + 1].id, lattice_list);

        h_target_lattice_2neighbour_temp[b * 6 + 2].id = lattice_list->getId(x4, y4, z4 - 1);
        h_target_lattice_2neighbour_temp[b * 6 + 2].type = getType(h_target_lattice_2neighbour_temp[b * 6 + 2].id, lattice_list);

        h_target_lattice_2neighbour_temp[b * 6 + 3].id = lattice_list->getId(x4, y4, z4 + 1);
        h_target_lattice_2neighbour_temp[b * 6 + 3].type = getType(h_target_lattice_2neighbour_temp[b * 6 + 3].id, lattice_list);

        h_target_lattice_2neighbour_temp[b * 6 + 4].id = lattice_list->getId(x4, y4 + 1, z4);
        h_target_lattice_2neighbour_temp[b * 6 + 4].type = getType(h_target_lattice_2neighbour_temp[b * 6 + 4].id, lattice_list);

        h_target_lattice_2neighbour_temp[b * 6 + 5].id = lattice_list->getId(x4 + 2, y4, z4);
        h_target_lattice_2neighbour_temp[b * 6 + 5].type = getType(h_target_lattice_2neighbour_temp[b * 6 + 5].id, lattice_list);

    }

    // end_time = MPI_Wtime() - start_time;
    // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", "  end_time 2 is : {} . \n", end_time);
    
    HANDLE_HIP(hipMemcpy(d_vacancy + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_lat_list + 8 * h_event_id, h_lat_list_temp, 8 * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    // HANDLE_HIP(hipMemcpy(d_lat_list_real + 8 * h_event_id, h_lat_list_temp_real, 8 * sizeof(unsigned long), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_2nn_neighbour + 6 * h_event_id, h_2nn_neighbour_temp, 6 * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_target_lattice_1neighbour + 64 * h_event_id, h_target_lattice_1neighbour_temp, 64 * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    HANDLE_HIP(hipMemcpy(d_target_lattice_2neighbour + 48 * h_event_id, h_target_lattice_2neighbour_temp, 48 * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
    
  }

  // dev_Vacancy *h_vacancy_temp;
  // h_vacancy_temp = (dev_Vacancy *)malloc(vac_count * sizeof(dev_Vacancy));
  // HANDLE_HIP(hipMemcpy(h_vacancy_temp, d_vacancy, vac_count * sizeof(dev_Vacancy), hipMemcpyDeviceToHost));
  // for(int i = 0; i < vac_count; i++){
  //   if(h_vacancy_temp[i].valid == 1){
  //     kiwi::logs::v(" ", " id is : {} .\n", h_vacancy_temp[i].id);
  //   }
  // }
  // free(h_vacancy_temp);
  sum_rates = 0;

  blockSize = 64;
  numBlocks = (vac_count + blockSize - 1) / blockSize;


  hipLaunchKernelGGL(updatennLatti, numBlocks, blockSize, 0, 0, 
                     d_selected,
                     d_lat_list,
                     d_2nn_neighbour, 
                     d_target_lattice_1neighbour,
                     d_target_lattice_2neighbour,
                     vac_count);
  CHECK_FOR_ERROR();
  // 调用GPU内核函数
  //start_time = MPI_Wtime();
  hipLaunchKernelGGL(calcRatesGPU, numBlocks, blockSize, 0, 0, 
                  d_vacancy,
                  d_lat_list,
                  d_2nn_neighbour, 
                  d_target_lattice_1neighbour,
                  d_target_lattice_2neighbour,
                  d_config,
                  d_sum_rates, 
                  vac_count);
  // hipDeviceSynchronize();
  // end_time = MPI_Wtime() - start_time;
  // total_time1 += end_time;
  CHECK_FOR_ERROR();

  HANDLE_HIP(hipMemcpy(&sum_rates, d_sum_rates, sizeof(double), hipMemcpyDeviceToHost));
  return sum_rates;
}

void sector_final(){
  using namespace gpuSimulate;

  HANDLE_HIP(hipFree(d_vacancy));
  HANDLE_HIP(hipFree(d_lat_list));
  // HANDLE_HIP(hipFree(d_lat_list_real));
  HANDLE_HIP(hipFree(d_2nn_neighbour));
  HANDLE_HIP(hipFree(d_target_lattice_1neighbour));
  HANDLE_HIP(hipFree(d_target_lattice_2neighbour));

  // free(h_vacancy_list);
}

void initialize_gpu(int process_rank){
    int gpu_count;
    HANDLE_HIP(hipGetDeviceCount(&gpu_count));
    HANDLE_HIP(hipSetDevice(process_rank % gpu_count));
    HANDLE_HIP(hipDeviceReset());
}

void gpu_prepare(const double v, const double T){
  using namespace gpuSimulate;

  //long h_size[2];
  double h_config[2];
  // h_size[0] = size_x;   // size_x
  // h_size[1] = size_y;   // size_y
  h_config[0] = v;   // v
  h_config[1] = T;   // T

  // HANDLE_HIP(hipMalloc((void**)&d_size, 2 * sizeof(long)));
  HANDLE_HIP(hipMalloc((void**)&d_config, 2 * sizeof(double)));

  //HANDLE_HIP(hipMemcpy(d_size, h_size, 2 * sizeof(long), hipMemcpyHostToDevice));
  HANDLE_HIP(hipMemcpy(d_config, h_config, 2 * sizeof(double), hipMemcpyHostToDevice));

  HANDLE_HIP(hipMalloc((void**)&d_sum_rates, sizeof(double)));
  // HANDLE_HIP(hipMalloc((void**)&d_event, sizeof(dev_event)));
  HANDLE_HIP(hipMalloc((void**)&d_selected, sizeof(dev_event)));
  HANDLE_HIP(hipMalloc((void**)&d_event_id, sizeof(int)));

  // HANDLE_HIP(hipMalloc((void**)&d_index, sizeof(int)));
  // HANDLE_HIP(hipMalloc((void**)&d_excepted_rate, sizeof(double)));
  //HANDLE_HIP(hipMalloc((void**)&d_selected_event_tag, sizeof(bool)));

  // HANDLE_HIP(hipMalloc((void**)&energy_temp1, 8 * sizeof(unsigned long)));
  // HANDLE_HIP(hipMalloc((void**)&energy_temp2, 8 * sizeof(LatticeTypes::lat_type)));

  // HANDLE_HIP(hipHostMalloc((void **)&h_lat_list_temp, 8 * sizeof(dev_nnLattice)));
  // HANDLE_HIP(hipHostMalloc((void **)&h_lat_list_temp_real, 8 * sizeof(unsigned long)));
  // HANDLE_HIP(hipHostMalloc((void **)&h_2nn_neighbour_temp, 6 * sizeof(dev_nnLattice)));
  // HANDLE_HIP(hipHostMalloc((void **)&h_target_lattice_1neighbour_temp, 64 * sizeof(dev_nnLattice)));
  // HANDLE_HIP(hipHostMalloc((void **)&h_target_lattice_2neighbour_temp, 48 * sizeof(dev_nnLattice)));

  h_lat_list_temp = (dev_nnLattice *)malloc(8 * sizeof(dev_nnLattice));
  // h_lat_list_temp_real = (unsigned long *)malloc(8 * sizeof(unsigned long));
  h_2nn_neighbour_temp = (dev_nnLattice *)malloc(6 * sizeof(dev_nnLattice));
  h_target_lattice_1neighbour_temp = (dev_nnLattice *)malloc(64 * sizeof(dev_nnLattice));
  h_target_lattice_2neighbour_temp = (dev_nnLattice *)malloc(48 * sizeof(dev_nnLattice));

  HANDLE_HIP(hipHostRegister(h_lat_list_temp, 8 * sizeof(dev_nnLattice), hipHostRegisterDefault));
  // HANDLE_HIP(hipHostRegister(h_lat_list_temp_real, 8 * sizeof(unsigned long), hipHostRegisterDefault));
  HANDLE_HIP(hipHostRegister(h_2nn_neighbour_temp, 6 * sizeof(dev_nnLattice), hipHostRegisterDefault));
  HANDLE_HIP(hipHostRegister(h_target_lattice_1neighbour_temp, 64 * sizeof(dev_nnLattice), hipHostRegisterDefault));
  HANDLE_HIP(hipHostRegister(h_target_lattice_2neighbour_temp, 48 * sizeof(dev_nnLattice), hipHostRegisterDefault));

}

void gpu_final(){

    using namespace gpuSimulate;
    // 释放GPU内存
    // HANDLE_HIP(hipFree(d_vac_rate));
    //HANDLE_HIP(hipFree(d_size));
    HANDLE_HIP(hipFree(d_config));
    // HANDLE_HIP(hipFree(d_region));
    // HANDLE_HIP(hipFree(d_box));
    HANDLE_HIP(hipFree(d_sum_rates));
    // HANDLE_HIP(hipFree(d_event));
    HANDLE_HIP(hipFree(d_selected));
    HANDLE_HIP(hipFree(d_event_id));

    // HANDLE_HIP(hipFree(d_index));
    // HANDLE_HIP(hipFree(d_excepted_rate));
    // HANDLE_HIP(hipHostFree(h_lat_list_temp));
    // HANDLE_HIP(hipHostFree(h_lat_list_temp_real));
    // HANDLE_HIP(hipHostFree(h_2nn_neighbour_temp));
    // HANDLE_HIP(hipHostFree(h_target_lattice_1neighbour_temp));
    // HANDLE_HIP(hipHostFree(h_target_lattice_2neighbour_temp));

    hipHostUnregister(h_lat_list_temp);
    // hipHostUnregister(h_lat_list_temp_real);
    hipHostUnregister(h_2nn_neighbour_temp);
    hipHostUnregister(h_target_lattice_1neighbour_temp);
    hipHostUnregister(h_target_lattice_2neighbour_temp);

    free(h_lat_list_temp);
    // free(h_lat_list_temp_real);
    free(h_2nn_neighbour_temp);
    free(h_target_lattice_1neighbour_temp);
    free(h_target_lattice_2neighbour_temp);

}

LatticeTypes getType(_type_lattice_id latti_id, LatticesList *lattice_list){
  auto it_vac = lattice_list->vac_hash.find(latti_id);
  auto it_cu = lattice_list->cu_hash.find(latti_id);
  auto it_mn = lattice_list->mn_hash.find(latti_id);
  auto it_ni = lattice_list->ni_hash.find(latti_id);
  if (it_vac != lattice_list->vac_hash.end()) {
    return LatticeTypes{LatticeTypes::V};
  }else if (it_cu != lattice_list->cu_hash.end()){
    return LatticeTypes{LatticeTypes::Cu};
  }else if (it_mn != lattice_list->mn_hash.end()){
    return LatticeTypes{LatticeTypes::Mn};
  }else if (it_ni != lattice_list->ni_hash.end()){
    return LatticeTypes{LatticeTypes::Ni};
  }else{
    return LatticeTypes{LatticeTypes::Fe};
  }
}

void exchange_lattices(LatticesList *lattice_list, _type_lattice_id from_id, _type_lattice_id to_id){
  auto it_vac = lattice_list->vac_hash.find(from_id);
  if (it_vac != lattice_list->vac_hash.end()) {
    lattice_list->vac_hash.erase(it_vac);
    lattice_list->vac_hash.emplace(std::make_pair(to_id, VacancyHash{}));
  }else{
    assert(false);
  }
  auto it_cu = lattice_list->cu_hash.find(to_id);
  auto it_mn = lattice_list->mn_hash.find(to_id);
  auto it_ni = lattice_list->ni_hash.find(to_id);
  if (it_cu != lattice_list->cu_hash.end()){
    lattice_list->cu_hash.erase(it_cu);
    lattice_list->cu_hash.emplace(from_id);
  }else if (it_mn != lattice_list->mn_hash.end()){
    lattice_list->mn_hash.erase(it_mn);
    lattice_list->mn_hash.emplace(from_id);
  }else if (it_ni != lattice_list->ni_hash.end()){
    lattice_list->ni_hash.erase(it_ni);
    lattice_list->ni_hash.emplace(from_id);
  }else{
    // Fe 不做操作
  }
}

void add_ghost(const comm::ColoredDomain *p_domain, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
               _type_lattice_id x, _type_lattice_id y, _type_lattice_id z, Lattice lat_to, const unsigned int sector_id){
      
  _type_lattice_id x_low, x_high, y_low, y_high, z_low, z_high;
  // 通过坐标判断处于哪一个 ghost 区中，一个子域有 7 个不同的ghost区
  // 7个exchange_ghost按照 不需要存储转发，需要存储转发1次(按照转发的顺序依次存储)，需要存储转发2次，依次存储
  switch (sector_id)
  {
    case 0:
      // 不需要存储转发 直接向下
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向前
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向下，再向前，最后向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向前
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向前，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 1:
      // 不需要存储转发 直接向下
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.x_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      //  需要存储转发1次 先向下，再向前
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.x_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向下，再向前，最后向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向前
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向前，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high); 
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 2:
      // 不需要存储转发 直接向下
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向后
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向下，再向后，最后向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向后
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向后，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 3:
      // 不需要存储转发 直接向下
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向后
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向下，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向下，再向后，最后向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low - p_domain->lattice_size_ghost[2]; // 0
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_low; // 2
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向后
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向后，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high); 
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_low; // 2
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_split_coord[2] + p_domain->lattice_size_ghost[2]; 
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 4:
      // 不需要存储转发 直接向上
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向前
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向上，再向前，最后向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向前
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向前，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 5:
      // 不需要存储转发 直接向上
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      //  需要存储转发1次 先向上，再向前
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向上，再向前，最后向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向前
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向前，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_low - p_domain->lattice_size_ghost[1]; // 0
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high); 
      y_low = p_domain->local_sub_box_lattice_region.y_low; // 2
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_split_coord[1] + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 6:
      // 不需要存储转发 直接向上
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向后
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向上，再向后，最后向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向后
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_split_coord[0] + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向后，再向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向左
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_low - p_domain->lattice_size_ghost[0]); // 0
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_low); // 2
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    case 7:
      // 不需要存储转发 直接向上
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[0].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向后
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[1].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向上，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[2].push_back(latti);
        break;
      }
      // 需要存储转发2次 先向上，再向后，最后向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_sub_box_lattice_region.z_high;
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high + p_domain->lattice_size_ghost[2];
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[3].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向后
      x_low = 2 * (p_domain->local_split_coord[0] - p_domain->lattice_size_ghost[0]);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[4].push_back(latti);
        break;
      }
      // 需要存储转发1次 先向后，再向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high);
      y_low = p_domain->local_sub_box_lattice_region.y_high;
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high + p_domain->lattice_size_ghost[1];
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[5].push_back(latti);
        break;
      }

      // 不需要存储转发 直接向右
      x_low = 2 * (p_domain->local_sub_box_lattice_region.x_high); 
      y_low = p_domain->local_split_coord[1] - p_domain->lattice_size_ghost[1];
      z_low = p_domain->local_split_coord[2] - p_domain->lattice_size_ghost[2]; 
      x_high = 2 * (p_domain->local_sub_box_lattice_region.x_high + p_domain->lattice_size_ghost[0]);
      y_high = p_domain->local_sub_box_lattice_region.y_high;
      z_high = p_domain->local_sub_box_lattice_region.z_high;
      if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
        ChangeLattice latti;
        latti.type = lat_to.type;
        latti.x = x;
        latti.y = y;
        latti.z = z;
        exchange_ghost[6].push_back(latti);
        break;
      }
      assert(false);
      break;
    default:
      assert(false);
      break;
  }
}