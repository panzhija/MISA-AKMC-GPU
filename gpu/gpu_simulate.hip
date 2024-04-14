#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include "gpu_simulate.h"
#include <map>
#include <iostream>

// Global Variables
namespace gpuSimulate {

    // gpu variable
    double *d_sum_rates;
    double *d_sum_rates_adv;
    double *d_sum_rates_aft;
    double *d_config;
    dev_event *d_event;
    _type_lattice_count *d_event_id;
    // int *d_selected; // 大小为 2 的数组，其第一个元素代表选中的事件中的空位是 h_vacancy_list 数组下的第 d_selected[0] 个元素，它向着它的 d_selected[1] 方向的 1nn跃迁.

    dev_Vacancy *d_vacancy_adv;
    dev_nnLattice *d_nnneighbour_adv;

    dev_Vacancy *d_vacancy_aft;
    dev_nnLattice *d_nnneighbour_aft;

    dev_Vacancy *d_vacancy_buf_adv;
    dev_nnLattice *d_nnneighbour_buf_adv;

    dev_Vacancy *d_vacancy_buf_aft;
    dev_nnLattice *d_nnneighbour_buf_aft;

    // host variable
    dev_Vacancy *h_vacancy_adv;
    dev_nnLattice *h_nnneighbour_adv;

    dev_Vacancy *h_vacancy_aft;
    dev_nnLattice *h_nnneighbour_aft;

    dev_Vacancy *h_vacancy_buf_adv;
    dev_nnLattice *h_nnneighbour_buf_adv;

    dev_Vacancy *h_vacancy_buf_aft;
    dev_nnLattice *h_nnneighbour_buf_aft;

    int blockSize;
    int numBlocks;
    _type_lattice_count vac_host_array_adv_size;
    _type_lattice_count vac_host_array_aft_size;
    _type_lattice_count vac_host_array_buf_adv_size;
    _type_lattice_count vac_host_array_buf_aft_size;

    // _type_lattice_id h_lat_list_id[8];

    double *sum_rates;
    double *sum_rates_adv;
    double *sum_rates_aft;
    _type_lattice_count vac_count_adv = 0;
    _type_lattice_count vac_count_aft = 0;
    _type_lattice_count vac_count_buf_adv = 0;
    _type_lattice_count vac_count_buf_aft = 0;
    _type_lattice_count vac_valid_count = 0;
    _type_lattice_count vac_valid_buf_count = 0;
    dev_event h_event;
    _type_lattice_count h_event_id;
    _type_lattice_count initial_size[4];

    _type_lattice_id to_x;
    _type_lattice_id to_y;
    _type_lattice_id to_z;
    int last_is_ghost;
    dev_Vacancy h_vacancy_temp;
    dev_nnLattice *h_nnneighbour_temp;

    hipStream_t stream[4];
    hipStream_t stream_buffer[4];
    LatticesList *lattice_list;
    comm::ColoredDomain *p_domain;
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

__global__ void calcRatesGPU(dev_event *d_event, dev_Vacancy *d_vacancy_adv, dev_Vacancy *d_vacancy_aft, 
                             dev_nnLattice *d_nnneighbour_adv, dev_nnLattice *d_nnneighbour_aft,
                             double *d_config, double* d_sum_rates, _type_lattice_count vac_count_adv, _type_lattice_count vac_total) {
    // 获取线程的全局ID
    _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
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
    if(tid < vac_count_adv && d_vacancy_adv[tid].valid != 0) {
      // 通过 d_event 更新由于上一次跃迁导致某些变化的邻居数据
      for(_type_lattice_count i = 0; i < NN_TOTAL; i++) {
        if(d_nnneighbour_adv[tid * NN_TOTAL + i].id == d_event->to_id) d_nnneighbour_adv[tid * NN_TOTAL + i].type._type = LatticeTypes::V;
        if(d_nnneighbour_adv[tid * NN_TOTAL + i].id == d_event->from_id) d_nnneighbour_adv[tid * NN_TOTAL + i].type._type = d_event->to_type._type;
      }
      d_vacancy_adv[tid].rates[0] = 0;
      d_vacancy_adv[tid].rates[1] = 0;
      d_vacancy_adv[tid].rates[2] = 0;
      d_vacancy_adv[tid].rates[3] = 0;
      d_vacancy_adv[tid].rates[4] = 0;
      d_vacancy_adv[tid].rates[5] = 0;
      d_vacancy_adv[tid].rates[6] = 0;
      d_vacancy_adv[tid].rates[7] = 0;

      double vac_energy = 0;
      int first_ir = 1;
      unsigned char atom_nei_status = 0;
      for (int b = 0; b < 8; b++) {
        if (d_nnneighbour_adv[tid * NN_TOTAL + b].type._type <= LatticeTypes::Mn && d_nnneighbour_adv[tid * NN_TOTAL + b].type._type != LatticeTypes::V) { // can trans
          atom_nei_status |= 1 << b;
        }
      }

      for (unsigned char b = 0; b < 8; b++) {
        if ((atom_nei_status >> b) & 1) { // the neighbour lattice is available
          double active_energy;
          switch (d_nnneighbour_adv[tid * NN_TOTAL + b].type._type) {
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
                  if (d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_1nn = d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type;
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

                  if (d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type;
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
                if (d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
                if (d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
                LatticeTypes::lat_type temp_type = d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type;
                if (b2 == b){
                  temp_type = LatticeTypes::V;
                }
                if (temp_type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = temp_type; 
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
                if (d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
                LatticeTypes::lat_type temp_type = d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
                if (b2 + b == 7) {
                  temp_type = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
                }
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
                if (d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;;
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
          d_vacancy_adv[tid].rates[b] = rate;
          // rate_total += rate;
          // return arrhenius(attempt_freq, temperature, active_energy);

          // _type_rate rate = d_vac_rate.rate(lattice, *lat_nei, ghost_atom, _1nn_offset);
          atomicAdd(&tmp_sh, rate);
        }
      }
      // d_vacancy[tid].rates[8] = rate_total;
      // atomicAdd(&tmp_sh, rate_total);
    }else if(vac_count_adv <= tid && tid < vac_total && d_vacancy_aft[tid - vac_count_adv].valid != 0) {

      // 通过 d_event 更新由于上一次跃迁导致某些变化的邻居数据
      for(_type_lattice_count i = 0; i < NN_TOTAL; i++){
        if(d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].id == d_event->to_id) d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].type._type = LatticeTypes::V;
        if(d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].id == d_event->from_id) d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].type._type = d_event->to_type._type;
      }

      d_vacancy_aft[tid - vac_count_adv].rates[0] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[1] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[2] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[3] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[4] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[5] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[6] = 0;
      d_vacancy_aft[tid - vac_count_adv].rates[7] = 0;

      double vac_energy = 0;
      int first_ir = 1;
      unsigned char atom_nei_status = 0;
      for (int b = 0; b < 8; b++) {
        if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type <= LatticeTypes::Mn && d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type != LatticeTypes::V) { // can trans
          atom_nei_status |= 1 << b;
        }
      }

      for (unsigned char b = 0; b < 8; b++) {
        if ((atom_nei_status >> b) & 1) { // the neighbour lattice is available
          double active_energy;
          switch (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type) {
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
                  if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_1nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b2].type._type;
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

                  if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_2nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 8 + b2].type._type;
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
                if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 14 + b * 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 14 + b * 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type;
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
                if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 78 + b * 6 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type;
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
                LatticeTypes::lat_type temp_type = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b2].type._type;
                if (b2 == b){
                  temp_type = LatticeTypes::V;
                }
                if (temp_type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = temp_type; 
                  LatticeTypes::lat_type temp_src = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type;
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
                if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type;
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
                LatticeTypes::lat_type temp_type = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 14 + b * 8 + b2].type._type;
                if (b2 + b == 7) {
                  temp_type = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + b].type._type;
                }
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
                if (d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + 78 + b * 6 + b2].type._type;;
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
          d_vacancy_aft[tid - vac_count_adv].rates[b] = rate;
          // rate_total += rate;
          // return arrhenius(attempt_freq, temperature, active_energy);

          // _type_rate rate = d_vac_rate.rate(lattice, *lat_nei, ghost_atom, _1nn_offset);
          atomicAdd(&tmp_sh, rate);
        }
      }
    }
  __syncthreads();
  if (threadIdx.x == 0) atomicAdd(d_sum_rates, tmp_sh);
}

__global__ void calcRatesGPUAdv(dev_Vacancy *d_vacancy_adv, dev_nnLattice *d_nnneighbour_adv,
                             double *d_config, double* d_sum_rates, _type_lattice_count vac_total) {
    // 获取线程的全局ID
    _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
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
    if(tid < vac_total) {
      d_vacancy_adv[tid].rates[0] = 0;
      d_vacancy_adv[tid].rates[1] = 0;
      d_vacancy_adv[tid].rates[2] = 0;
      d_vacancy_adv[tid].rates[3] = 0;
      d_vacancy_adv[tid].rates[4] = 0;
      d_vacancy_adv[tid].rates[5] = 0;
      d_vacancy_adv[tid].rates[6] = 0;
      d_vacancy_adv[tid].rates[7] = 0;

      double vac_energy = 0;
      int first_ir = 1;
      unsigned char atom_nei_status = 0;
      for (int b = 0; b < 8; b++) {
        if (d_nnneighbour_adv[tid * NN_TOTAL + b].type._type <= LatticeTypes::Mn && d_nnneighbour_adv[tid * NN_TOTAL + b].type._type != LatticeTypes::V) { // can trans
          atom_nei_status |= 1 << b;
        }
      }

      for (unsigned char b = 0; b < 8; b++) {
        if ((atom_nei_status >> b) & 1) { // the neighbour lattice is available
          double active_energy;
          switch (d_nnneighbour_adv[tid * NN_TOTAL + b].type._type) {
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
                  if (d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_1nn = d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type;
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

                  if (d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                    // it is vacancy or single atom.
                    LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type;
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

            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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

            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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

// ______________________________________________________________________________________________________________________________________



          // calculate system energy after transition.

          double e_after = 0;
          {
            // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
            double e_src, e_des;
            unsigned char _1nn_status = 0xFF;
                
            double energy = 0.0;

            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                LatticeTypes::lat_type temp_type = d_nnneighbour_adv[tid * NN_TOTAL + b2].type._type;
                if (b2 == b){
                  temp_type = LatticeTypes::V;
                }
                if (temp_type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = temp_type; 
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
                if (d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
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
            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                LatticeTypes::lat_type temp_type = d_nnneighbour_adv[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
                if (b2 + b == 7) {
                  temp_type = d_nnneighbour_adv[tid * NN_TOTAL + b].type._type;
                }
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
            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Mn) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = d_nnneighbour_adv[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;;
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
          d_vacancy_adv[tid].rates[b] = rate;
          atomicAdd(&tmp_sh, rate);
        }
      }
    }
  __syncthreads();
  if (threadIdx.x == 0) atomicAdd(d_sum_rates, tmp_sh);
}

__global__ void selectEventGPU(dev_Vacancy *d_vac_adv, dev_Vacancy *d_vac_aft, dev_nnLattice *d_nn_adv, 
                               dev_nnLattice *d_nn_aft, dev_event *d_event, _type_lattice_count *d_event_id, _type_lattice_count vac_adv, 
                               _type_lattice_count vac_total, double excepted_rand) {

  // 获取线程的全局ID
  _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid == 0){
    __shared__ double excepted_rate;
    __shared__ double rate_accumulator;
    excepted_rate = excepted_rand;
    rate_accumulator = 0;

    int tag = 1;
    for(_type_lattice_count i = 0; i < vac_adv; i++){
      if(d_vac_adv[i].valid != 0) {
        if(tag == 1) {
          for (int rate_index = 0; rate_index < 8; rate_index++) {
            rate_accumulator += d_vac_adv[i].rates[rate_index];
            if(rate_accumulator > excepted_rate) {
              if(tag == 1) {
                d_event->from_id = d_vac_adv[i].id;
                d_event->to_id = d_nn_adv[i * NN_TOTAL + rate_index].id;
                d_event->to_type._type = d_nn_adv[i * NN_TOTAL + rate_index].type._type;
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

    for(_type_lattice_count i = vac_adv; i < vac_total; i++){
      if(d_vac_aft[i - vac_adv].valid != 0) {
        if(tag == 1) {
          for (int rate_index = 0; rate_index < 8; rate_index++) {
            rate_accumulator += d_vac_aft[i - vac_adv].rates[rate_index];
            if(rate_accumulator > excepted_rate) {
              if(tag == 1) {
                d_event->from_id = d_vac_aft[i - vac_adv].id;
                d_event->to_id = d_nn_aft[(i - vac_adv) * NN_TOTAL + rate_index].id;
                d_event->to_type._type = d_nn_aft[(i - vac_adv) * NN_TOTAL + rate_index].type._type;
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


// __global__ void updatennLatti(dev_event *d_event, dev_nnLattice *d_nnneighbour_adv,
//                               dev_nnLattice *d_nnneighbour_aft, _type_lattice_count vac_count_adv, _type_lattice_count vac_total) {
//   _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
//   if(tid < vac_count_adv) {
//     for(_type_lattice_count i = 0; i < NN_TOTAL; i++){
//       if(d_nnneighbour_adv[tid * NN_TOTAL + i].id == d_event->to_id) d_nnneighbour_adv[tid * NN_TOTAL + i].type._type = LatticeTypes::V;
//       if(d_nnneighbour_adv[tid * NN_TOTAL + i].id == d_event->from_id) d_nnneighbour_adv[tid * NN_TOTAL + i].type._type = d_event->to_type._type;
//     }
//   }else if(vac_count_adv <= tid && tid < vac_total) {
//     for(_type_lattice_count i = 0; i < NN_TOTAL; i++){
//       if(d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].id == d_event->to_id) d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].type._type = LatticeTypes::V;
//       if(d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].id == d_event->from_id) d_nnneighbour_aft[(tid - vac_count_adv) * NN_TOTAL + i].type._type = d_event->to_type._type;
//     }
//   }
// }

void firstSectLocalToGpu(_type_lattice_id *vac_idArray_adv, _type_lattice_count vac_adv) {
  using namespace gpuSimulate;

  vac_count_adv = vac_adv;
  vac_valid_count += vac_adv;
  
  int i = 0;
  if(vac_count_adv > 0) {
    _type_lattice_count index = 0;
    for(i = 0; i < 4; i++){
      initial_size[i] = vac_count_adv / 4 + (i < (vac_count_adv % 4) ? 1 : 0);
      if(initial_size[i] > 0) {
        initInfo(initial_size[i], vac_idArray_adv + index, h_vacancy_adv + index, h_nnneighbour_adv + index * NN_TOTAL);
        HANDLE_HIP(hipMemcpyAsync(d_vacancy_adv + index, h_vacancy_adv + index, initial_size[i] * sizeof(dev_Vacancy), hipMemcpyHostToDevice, stream[i]));
        HANDLE_HIP(hipMemcpyAsync(d_nnneighbour_adv + index * NN_TOTAL, h_nnneighbour_adv + index * NN_TOTAL, NN_TOTAL * initial_size[i] * sizeof(dev_nnLattice), hipMemcpyHostToDevice, stream[i]));
        index += initial_size[i];
      }else{
        break;
      }
    }
    // initInfo(vac_count_adv, vac_idArray_adv, h_vacancy_adv, h_nnneighbour_adv);

    // HANDLE_HIP(hipMemcpy(d_vacancy_adv, h_vacancy_adv, vac_count_adv * sizeof(dev_Vacancy), hipMemcpyHostToDevice));
    // HANDLE_HIP(hipMemcpy(d_nnneighbour_adv, h_nnneighbour_adv, NN_TOTAL * vac_count_adv * sizeof(dev_nnLattice), hipMemcpyHostToDevice));

    blockSize = 64;
    numBlocks = (vac_count_adv + blockSize - 1) / blockSize;

    // 调用GPU内核函数
    hipLaunchKernelGGL(calcRatesGPUAdv, numBlocks, blockSize, 0, stream[i - 1], 
                       d_vacancy_adv,
                       d_nnneighbour_adv, 
                       d_config, 
                       d_sum_rates_adv, 
                       vac_count_adv);
    // free(h_vacancy_adv);
    // free(h_nnneighbour_adv);
  }
  free(vac_idArray_adv);
  // kiwi::logs::v(" ", " first vac_adv count GPU is : {}.\n", vac_adv);
}

double first_calculate_GPU(_type_lattice_id *vac_idArray_aft, _type_lattice_count vac_aft, _type_lattice_id *vac_idArray_buf_adv, _type_lattice_count vac_adv, int sect){
  
  using namespace gpuSimulate;
    
    // size_t free_mem, total_mem;
    // HANDLE_HIP(hipMemGetInfo(&free_mem, &total_mem));
    // float free_memory = free_mem / (1024.0 * 1024.0);
    // float total_memory = total_mem / (1024.0 * 1024.0); 
    // kiwi::logs::v(" ", "Device_mem_before: free = {} total = {}.\n", free_memory, total_memory);
    // for(int i = 0; i < vac_total; i++){
    //   kiwi::logs::v(" ", " id is : {} x is : {} y is : {} z is : {}.\n", vac_total[i]);
    // }
    
    *sum_rates = 0;
    *sum_rates_adv = 0;
    *sum_rates_aft = 0;
    int aft_i = 0;
    if(sect % 2 == 0){
      // kiwi::logs::v(" ", " vac_aft count GPU is : {}.\n", vac_aft);
      // kiwi::logs::v(" ", " next_vac_buf_adv count GPU is : {}.\n", vac_adv);
      vac_count_aft = vac_aft;
      vac_valid_count += vac_aft;
      vac_count_buf_adv = vac_adv;
      vac_valid_buf_count += vac_adv;

      // 重新分配锁页内存和 device 内存
      if(vac_count_aft > vac_host_array_aft_size){
        vac_host_array_aft_size = vac_count_aft * RE_SCALE_SIZE;
        hipHostFree(h_vacancy_aft);
        hipHostFree(h_nnneighbour_aft);
        HANDLE_HIP(hipFree(d_vacancy_aft));
        HANDLE_HIP(hipFree(d_nnneighbour_aft));
        hipHostMalloc((void**)&h_vacancy_aft, vac_host_array_aft_size * sizeof(dev_Vacancy));
        hipHostMalloc((void**)&h_nnneighbour_aft, vac_host_array_aft_size * NN_TOTAL * sizeof(dev_nnLattice));
        HANDLE_HIP(hipMalloc((void**)&d_vacancy_aft, vac_host_array_aft_size * sizeof(dev_Vacancy)));
        HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_aft, vac_host_array_aft_size * NN_TOTAL * sizeof(dev_nnLattice)));
      }
      if(vac_count_aft > 0) {
        _type_lattice_count index = 0;
        for(aft_i = 0; aft_i < 4; aft_i++){
          initial_size[aft_i] = vac_count_aft / 4 + (aft_i < (vac_count_aft % 4) ? 1 : 0);
          if(initial_size[aft_i] > 0) {
            initInfo(initial_size[aft_i], vac_idArray_aft + index, h_vacancy_aft + index, h_nnneighbour_aft + index * NN_TOTAL);
            HANDLE_HIP(hipMemcpyAsync(d_vacancy_aft + index, h_vacancy_aft + index, initial_size[aft_i] * sizeof(dev_Vacancy), hipMemcpyHostToDevice, stream[aft_i]));
            HANDLE_HIP(hipMemcpyAsync(d_nnneighbour_aft + index * NN_TOTAL, h_nnneighbour_aft + index * NN_TOTAL, NN_TOTAL * initial_size[aft_i] * sizeof(dev_nnLattice), hipMemcpyHostToDevice, stream[aft_i]));
            index += initial_size[aft_i];
          }else{
            break;
          }
        }
        blockSize = 64;
        numBlocks = (vac_count_aft + blockSize - 1) / blockSize;

        // 调用GPU内核函数
        hipLaunchKernelGGL(calcRatesGPUAdv, numBlocks, blockSize, 0, stream[aft_i - 1], 
                           d_vacancy_aft,
                           d_nnneighbour_aft, 
                           d_config, 
                           d_sum_rates_aft, 
                           vac_count_aft);
      }
      //if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " sect is : {} calculate vac count is : {} .\n", sect, vac_count_adv + vac_count_aft);

      // 重新分配锁页内存和 device 内存
      if(vac_count_buf_adv > vac_host_array_buf_adv_size){
        vac_host_array_buf_adv_size = vac_count_buf_adv * RE_SCALE_SIZE;
        hipHostFree(h_vacancy_buf_adv);
        hipHostFree(h_nnneighbour_buf_adv);
        HANDLE_HIP(hipFree(d_vacancy_buf_adv));
        HANDLE_HIP(hipFree(d_nnneighbour_buf_adv));
        hipHostMalloc((void**)&h_vacancy_buf_adv, vac_host_array_buf_adv_size * sizeof(dev_Vacancy));
        hipHostMalloc((void**)&h_nnneighbour_buf_adv, vac_host_array_buf_adv_size * NN_TOTAL * sizeof(dev_nnLattice));
        HANDLE_HIP(hipMalloc((void**)&d_vacancy_buf_adv, vac_host_array_buf_adv_size * sizeof(dev_Vacancy)));
        HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_buf_adv, vac_host_array_buf_adv_size * NN_TOTAL * sizeof(dev_nnLattice)));
      }
      if(vac_count_buf_adv > 0){
        _type_lattice_count index = 0;
        for(int i = 0; i < 4; i++){
          initial_size[i] = vac_count_buf_adv / 4 + (i < (vac_count_buf_adv % 4) ? 1 : 0);
          if(initial_size[i] > 0) {
            initInfo(initial_size[i], vac_idArray_buf_adv + index, h_vacancy_buf_adv + index, h_nnneighbour_buf_adv + index * NN_TOTAL);
            HANDLE_HIP(hipMemcpyAsync(d_vacancy_buf_adv + index, h_vacancy_buf_adv + index, initial_size[i] * sizeof(dev_Vacancy), hipMemcpyHostToDevice, stream_buffer[i]));
            HANDLE_HIP(hipMemcpyAsync(d_nnneighbour_buf_adv + index * NN_TOTAL, h_nnneighbour_buf_adv + index * NN_TOTAL, NN_TOTAL * initial_size[i] * sizeof(dev_nnLattice), hipMemcpyHostToDevice, stream_buffer[i]));
            index += initial_size[i];
          }else{
            break;
          }
        }
      }
      if(vac_count_adv > 0) HANDLE_HIP(hipMemcpyAsync(sum_rates_adv, d_sum_rates_adv, sizeof(double), hipMemcpyDeviceToHost, stream[aft_i - 1]));
      if(vac_count_aft > 0) HANDLE_HIP(hipMemcpyAsync(sum_rates_aft, d_sum_rates_aft, sizeof(double), hipMemcpyDeviceToHost, stream[aft_i - 1]));
    }else{
      // kiwi::logs::v(" ", " vac_buf_aft count GPU is : {}.\n", vac_aft);
      // kiwi::logs::v(" ", " next_vac_adv count GPU is : {}.\n", vac_adv);
      vac_count_buf_aft = vac_aft;
      vac_valid_buf_count += vac_aft;
      vac_count_adv = vac_adv;
      vac_valid_count += vac_adv;

      if(vac_count_buf_aft > vac_host_array_buf_aft_size){
        vac_host_array_buf_aft_size = vac_count_buf_aft * RE_SCALE_SIZE;
        hipHostFree(h_vacancy_buf_aft);
        hipHostFree(h_nnneighbour_buf_aft);
        HANDLE_HIP(hipFree(d_vacancy_buf_aft));
        HANDLE_HIP(hipFree(d_nnneighbour_buf_aft));
        hipHostMalloc((void**)&h_vacancy_buf_aft, vac_host_array_buf_aft_size * sizeof(dev_Vacancy));
        hipHostMalloc((void**)&h_nnneighbour_buf_aft, vac_host_array_buf_aft_size * NN_TOTAL * sizeof(dev_nnLattice));
        HANDLE_HIP(hipMalloc((void**)&d_vacancy_buf_aft, vac_host_array_buf_aft_size * sizeof(dev_Vacancy)));
        HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_buf_aft, vac_host_array_buf_aft_size * NN_TOTAL * sizeof(dev_nnLattice)));
      }
      if(vac_count_buf_aft > 0) {
        _type_lattice_count index = 0;
        for(aft_i = 0; aft_i < 4; aft_i++){
          initial_size[aft_i] = vac_count_buf_aft / 4 + (aft_i < (vac_count_buf_aft % 4) ? 1 : 0);
          if(initial_size[aft_i] > 0) {
            initInfo(initial_size[aft_i], vac_idArray_aft + index, h_vacancy_buf_aft + index, h_nnneighbour_buf_aft + index * NN_TOTAL);
            HANDLE_HIP(hipMemcpyAsync(d_vacancy_buf_aft + index, h_vacancy_buf_aft + index, initial_size[aft_i] * sizeof(dev_Vacancy), hipMemcpyHostToDevice, stream_buffer[aft_i]));
            HANDLE_HIP(hipMemcpyAsync(d_nnneighbour_buf_aft + index * NN_TOTAL, h_nnneighbour_buf_aft + index * NN_TOTAL, NN_TOTAL * initial_size[aft_i] * sizeof(dev_nnLattice), hipMemcpyHostToDevice, stream_buffer[aft_i]));
            index += initial_size[aft_i];
          }else{
            break;
          }
        }
        blockSize = 64;
        numBlocks = (vac_count_buf_aft + blockSize - 1) / blockSize;

        // 调用GPU内核函数
        hipLaunchKernelGGL(calcRatesGPUAdv, numBlocks, blockSize, 0, stream_buffer[aft_i - 1], 
                           d_vacancy_buf_aft,
                           d_nnneighbour_buf_aft,
                           d_config,
                           d_sum_rates_aft, 
                           vac_count_buf_aft);
      }
      //if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " sect is : {} calculate vac count is : {} .\n", sect, vac_count_buf_adv + vac_count_buf_aft);
      if(vac_count_adv > vac_host_array_adv_size){
        vac_host_array_adv_size = vac_count_adv * RE_SCALE_SIZE;
        hipHostFree(h_vacancy_adv);
        hipHostFree(h_nnneighbour_adv);
        HANDLE_HIP(hipFree(d_vacancy_adv));
        HANDLE_HIP(hipFree(d_nnneighbour_adv));
        hipHostMalloc((void**)&h_vacancy_adv, vac_host_array_adv_size * sizeof(dev_Vacancy));
        hipHostMalloc((void**)&h_nnneighbour_adv, vac_host_array_adv_size * NN_TOTAL * sizeof(dev_nnLattice));
        HANDLE_HIP(hipMalloc((void**)&d_vacancy_adv, vac_host_array_adv_size * sizeof(dev_Vacancy)));
        HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_adv, vac_host_array_adv_size * NN_TOTAL * sizeof(dev_nnLattice)));
      }
      if(vac_count_adv > 0){
        _type_lattice_count index = 0;
        for(int i = 0; i < 4; i++){
          initial_size[i] = vac_count_adv / 4 + (i < (vac_count_adv % 4) ? 1 : 0);
          if(initial_size[i] > 0) {
            initInfo(initial_size[i], vac_idArray_buf_adv + index, h_vacancy_adv + index, h_nnneighbour_adv + index * NN_TOTAL);
            HANDLE_HIP(hipMemcpyAsync(d_vacancy_adv + index, h_vacancy_adv + index, initial_size[i] * sizeof(dev_Vacancy), hipMemcpyHostToDevice, stream[i]));
            HANDLE_HIP(hipMemcpyAsync(d_nnneighbour_adv + index * NN_TOTAL, h_nnneighbour_adv + index * NN_TOTAL, NN_TOTAL * initial_size[i] * sizeof(dev_nnLattice), hipMemcpyHostToDevice, stream[i]));
            index += initial_size[i];
          }else{
            break;
          }
        }
      }
      if(vac_count_buf_adv > 0) HANDLE_HIP(hipMemcpyAsync(sum_rates_adv, d_sum_rates_adv, sizeof(double), hipMemcpyDeviceToHost, stream_buffer[aft_i - 1]));
      if(vac_count_buf_aft > 0) HANDLE_HIP(hipMemcpyAsync(sum_rates_aft, d_sum_rates_aft, sizeof(double), hipMemcpyDeviceToHost, stream_buffer[aft_i - 1]));
    }
    // HANDLE_HIP(hipMemGetInfo(&free_mem, &total_mem));
    // free_memory = free_mem / (1024.0 * 1024.0);
    // total_memory = total_mem / (1024.0 * 1024.0);
    // kiwi::logs::v(" ", "Device_mem_after: free = {} total = {}.\n", free_memory, total_memory);
    CHECK_FOR_ERROR();
    free(vac_idArray_aft);
    free(vac_idArray_buf_adv);

    hipDeviceSynchronize();
    *sum_rates = *sum_rates_adv + *sum_rates_aft;

    return *sum_rates;
}

void selectAndPerformEventGPU(double excepted_rand, int rank, _type_lattice_count step, int sect,
                              std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, const unsigned int sector_id,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z){
  using namespace gpuSimulate;

  blockSize = 1;
  numBlocks = 1;

  if(sect % 2 == 0) {
    hipLaunchKernelGGL(selectEventGPU, numBlocks, blockSize, 0, 0, 
                       d_vacancy_adv,
                       d_vacancy_aft,
                       d_nnneighbour_adv,
                       d_nnneighbour_aft,
                       d_event,
                       d_event_id,
                       vac_count_adv,
                       vac_count_adv + vac_count_aft,
                       *sum_rates * excepted_rand);
  }else{
    hipLaunchKernelGGL(selectEventGPU, numBlocks, blockSize, 0, 0, 
                       d_vacancy_buf_adv,
                       d_vacancy_buf_aft,
                       d_nnneighbour_buf_adv,
                       d_nnneighbour_buf_aft,
                       d_event,
                       d_event_id,
                       vac_count_buf_adv,
                       vac_count_buf_adv + vac_count_buf_aft,
                       *sum_rates * excepted_rand);

  }
  CHECK_FOR_ERROR();
  //kiwi::logs::v("", "excepted_rate is : {} .\n", sum_rates * excepted_rand);

  HANDLE_HIP(hipMemcpy(&h_event, d_event, sizeof(dev_event), hipMemcpyDeviceToHost));

  last_is_ghost = 0;

  _type_lattice_id from_x = h_event.from_id % lattice_list->meta.size_x;
  _type_lattice_id from_y = (h_event.from_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  _type_lattice_id from_z = h_event.from_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

  to_x = h_event.to_id % lattice_list->meta.size_x;
  to_y = (h_event.to_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  to_z = h_event.to_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

  exchange_lattices(h_event);
  if(lattice_list->meta.isGhostLat(h_event.to_id)){
    last_is_ghost = 1;
    Lattice lat_to;
    lat_to.id = h_event.to_id;
    lat_to.type._type = LatticeTypes::V;
    add_ghost(exchange_ghost, to_x, to_y, to_z, lat_to, sector_id);
  }
  if (lattice_list->meta.isSurfaceLat(h_event.from_id)) add_surface(h_event.from_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
  if (lattice_list->meta.isSurfaceLat(h_event.to_id)) add_surface(h_event.to_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
  // add_surface(h_event.from_id, h_event.to_id, next_sector_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
  // kiwi::logs::v(" ", " from_id is : {} to_id is : {} .\n", h_event.from_id, h_event.to_id);
}

double calculate_GPU(const comm::Region<comm::_type_lattice_size> region, int sect){
  using namespace gpuSimulate;

  HANDLE_HIP(hipMemcpy(&h_event_id, d_event_id, sizeof(_type_lattice_count), hipMemcpyDeviceToHost));
  // 判断转移后的位置是否位于阴影区 或者 是否还位于当前子域内。
  if (last_is_ghost == 1 || !(2 * region.x_low <= to_x && to_x < 2 * region.x_high && region.y_low <= to_y && to_y < region.y_high && region.z_low <= to_z && to_z < region.z_high)){
    // 转移后的位置位于阴影区 或者 位于其他子域内。
    if(sect % 2 == 0) {
      vac_valid_count--;
      if(vac_valid_count > 0){
        // 当前子域还有空位
        h_vacancy_temp.valid = 0;
        if(h_event_id < vac_count_adv){
          HANDLE_HIP(hipMemcpy(d_vacancy_adv + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        }else{
          HANDLE_HIP(hipMemcpy(d_vacancy_aft + h_event_id - vac_count_adv, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        }
      }else{
        // 当前子域无空位
        return 0;
      }
    }else{
      vac_valid_buf_count--;
      if(vac_valid_buf_count > 0){
        // 当前子域还有空位
        h_vacancy_temp.valid = 0;
        if(h_event_id < vac_count_buf_adv){
          HANDLE_HIP(hipMemcpy(d_vacancy_buf_adv + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        }else{
          HANDLE_HIP(hipMemcpy(d_vacancy_buf_aft + h_event_id - vac_count_buf_adv, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        }
      }else{
        // 当前子域无空位
        return 0;
      }
    }
    
  
  }else{
    // 转移后的位置没有位于阴影区 并且 还在当前子域内。
    // 更新相关数据。
    h_vacancy_temp.id = h_event.to_id;
    h_vacancy_temp.valid = 1;

    if(to_x % 2 == 0){
      
      h_nnneighbour_temp[0].id = lattice_list->getId(to_x - 1, to_y - 1, to_z - 1);
      h_nnneighbour_temp[0].type._type = getType(h_nnneighbour_temp[0].id);

      h_nnneighbour_temp[1].id = lattice_list->getId(to_x - 1, to_y - 1, to_z);
      h_nnneighbour_temp[1].type._type = getType(h_nnneighbour_temp[1].id);

      h_nnneighbour_temp[2].id = lattice_list->getId(to_x - 1, to_y, to_z - 1);
      h_nnneighbour_temp[2].type._type = getType(h_nnneighbour_temp[2].id);

      h_nnneighbour_temp[3].id = lattice_list->getId(to_x - 1, to_y, to_z);
      h_nnneighbour_temp[3].type._type = getType(h_nnneighbour_temp[3].id);

      h_nnneighbour_temp[4].id = lattice_list->getId(to_x + 1, to_y - 1, to_z - 1);
      h_nnneighbour_temp[4].type._type = getType(h_nnneighbour_temp[4].id);

      h_nnneighbour_temp[5].id = lattice_list->getId(to_x + 1, to_y - 1, to_z);
      h_nnneighbour_temp[5].type._type = getType(h_nnneighbour_temp[5].id);

      h_nnneighbour_temp[6].id = lattice_list->getId(to_x + 1, to_y, to_z - 1);
      h_nnneighbour_temp[6].type._type = getType(h_nnneighbour_temp[6].id);

      h_nnneighbour_temp[7].id = lattice_list->getId(to_x + 1, to_y, to_z);
      h_nnneighbour_temp[7].type._type = getType(h_nnneighbour_temp[7].id);

    }else{

      h_nnneighbour_temp[0].id = lattice_list->getId(to_x - 1, to_y, to_z);
      h_nnneighbour_temp[0].type._type = getType(h_nnneighbour_temp[0].id);

      h_nnneighbour_temp[1].id = lattice_list->getId(to_x - 1, to_y, to_z + 1);
      h_nnneighbour_temp[1].type._type = getType(h_nnneighbour_temp[1].id);

      h_nnneighbour_temp[2].id = lattice_list->getId(to_x - 1, to_y + 1, to_z);
      h_nnneighbour_temp[2].type._type = getType(h_nnneighbour_temp[2].id);

      h_nnneighbour_temp[3].id = lattice_list->getId(to_x - 1, to_y + 1, to_z + 1);
      h_nnneighbour_temp[3].type._type = getType(h_nnneighbour_temp[3].id);

      h_nnneighbour_temp[4].id = lattice_list->getId(to_x + 1, to_y, to_z);
      h_nnneighbour_temp[4].type._type = getType(h_nnneighbour_temp[4].id);

      h_nnneighbour_temp[5].id = lattice_list->getId(to_x + 1, to_y, to_z + 1);
      h_nnneighbour_temp[5].type._type = getType(h_nnneighbour_temp[5].id);

      h_nnneighbour_temp[6].id = lattice_list->getId(to_x + 1, to_y + 1, to_z);
      h_nnneighbour_temp[6].type._type = getType(h_nnneighbour_temp[6].id);

      h_nnneighbour_temp[7].id = lattice_list->getId(to_x + 1, to_y + 1, to_z + 1);
      h_nnneighbour_temp[7].type._type = getType(h_nnneighbour_temp[7].id);

    }

    h_nnneighbour_temp[8].id = lattice_list->getId(to_x - 2, to_y, to_z);
    h_nnneighbour_temp[8].type._type = getType(h_nnneighbour_temp[8].id);

    h_nnneighbour_temp[9].id = lattice_list->getId(to_x, to_y - 1, to_z);
    h_nnneighbour_temp[9].type._type = getType(h_nnneighbour_temp[9].id);

    h_nnneighbour_temp[10].id = lattice_list->getId(to_x, to_y, to_z - 1);
    h_nnneighbour_temp[10].type._type = getType(h_nnneighbour_temp[10].id);

    h_nnneighbour_temp[11].id = lattice_list->getId(to_x, to_y, to_z + 1);
    h_nnneighbour_temp[11].type._type = getType(h_nnneighbour_temp[11].id);

    h_nnneighbour_temp[12].id = lattice_list->getId(to_x, to_y + 1, to_z);
    h_nnneighbour_temp[12].type._type = getType(h_nnneighbour_temp[12].id);

    h_nnneighbour_temp[13].id = lattice_list->getId(to_x + 2, to_y, to_z);
    h_nnneighbour_temp[13].type._type = getType(h_nnneighbour_temp[13].id);
        
    for(int b = 0; b < 8; b++){
      _type_lattice_id x4 = h_nnneighbour_temp[b].id % lattice_list->meta.size_x;
      _type_lattice_id y4 = (h_nnneighbour_temp[b].id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
      _type_lattice_id z4 = h_nnneighbour_temp[b].id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

      if(x4 % 2 == 0){
        
        h_nnneighbour_temp[14 + b * 8].id = lattice_list->getId(x4 - 1, y4 - 1, z4 - 1);
        h_nnneighbour_temp[14 + b * 8].type._type = getType(h_nnneighbour_temp[14 + b * 8].id);

        h_nnneighbour_temp[14 + b * 8 + 1].id = lattice_list->getId(x4 - 1, y4 - 1, z4);
        h_nnneighbour_temp[14 + b * 8 + 1].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 1].id);

        h_nnneighbour_temp[14 + b * 8 + 2].id = lattice_list->getId(x4 - 1, y4, z4 - 1);
        h_nnneighbour_temp[14 + b * 8 + 2].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 2].id);

        h_nnneighbour_temp[14 + b * 8 + 3].id = lattice_list->getId(x4 - 1, y4, z4);
        h_nnneighbour_temp[14 + b * 8 + 3].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 3].id);

        h_nnneighbour_temp[14 + b * 8 + 4].id = lattice_list->getId(x4 + 1, y4 - 1, z4 - 1);
        h_nnneighbour_temp[14 + b * 8 + 4].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 4].id);

        h_nnneighbour_temp[14 + b * 8 + 5].id = lattice_list->getId(x4 + 1, y4 - 1, z4);
        h_nnneighbour_temp[14 + b * 8 + 5].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 5].id);

        h_nnneighbour_temp[14 + b * 8 + 6].id = lattice_list->getId(x4 + 1, y4, z4 - 1);
        h_nnneighbour_temp[14 + b * 8 + 6].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 6].id);

        h_nnneighbour_temp[14 + b * 8 + 7].id = lattice_list->getId(x4 + 1, y4, z4);
        h_nnneighbour_temp[14 + b * 8 + 7].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 7].id);
          
      }else{

        h_nnneighbour_temp[14 + b * 8].id = lattice_list->getId(x4 - 1, y4, z4);
        h_nnneighbour_temp[14 + b * 8].type._type = getType(h_nnneighbour_temp[14 + b * 8].id);

        h_nnneighbour_temp[14 + b * 8 + 1].id = lattice_list->getId(x4 - 1, y4, z4 + 1);
        h_nnneighbour_temp[14 + b * 8 + 1].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 1].id);

        h_nnneighbour_temp[14 + b * 8 + 2].id = lattice_list->getId(x4 - 1, y4 + 1, z4);
        h_nnneighbour_temp[14 + b * 8 + 2].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 2].id);

        h_nnneighbour_temp[14 + b * 8 + 3].id = lattice_list->getId(x4 - 1, y4 + 1, z4 + 1);
        h_nnneighbour_temp[14 + b * 8 + 3].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 3].id);

        h_nnneighbour_temp[14 + b * 8 + 4].id = lattice_list->getId(x4 + 1, y4, z4);
        h_nnneighbour_temp[14 + b * 8 + 4].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 4].id);

        h_nnneighbour_temp[14 + b * 8 + 5].id = lattice_list->getId(x4 + 1, y4, z4 + 1);
        h_nnneighbour_temp[14 + b * 8 + 5].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 5].id);

        h_nnneighbour_temp[14 + b * 8 + 6].id = lattice_list->getId(x4 + 1, y4 + 1, z4);
        h_nnneighbour_temp[14 + b * 8 + 6].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 6].id);

        h_nnneighbour_temp[14 + b * 8 + 7].id = lattice_list->getId(x4 + 1, y4 + 1, z4 + 1);
        h_nnneighbour_temp[14 + b * 8 + 7].type._type = getType(h_nnneighbour_temp[14 + b * 8 + 7].id);

      }

      h_nnneighbour_temp[78 + b * 6].id = lattice_list->getId(x4 - 2, y4, z4);
      h_nnneighbour_temp[78 + b * 6].type._type = getType(h_nnneighbour_temp[78 + b * 6].id);

      h_nnneighbour_temp[78 + b * 6 + 1].id = lattice_list->getId(x4, y4 - 1, z4);
      h_nnneighbour_temp[78 + b * 6 + 1].type._type = getType(h_nnneighbour_temp[78 + b * 6 + 1].id);

      h_nnneighbour_temp[78 + b * 6 + 2].id = lattice_list->getId(x4, y4, z4 - 1);
      h_nnneighbour_temp[78 + b * 6 + 2].type._type = getType(h_nnneighbour_temp[78 + b * 6 + 2].id);

      h_nnneighbour_temp[78 + b * 6 + 3].id = lattice_list->getId(x4, y4, z4 + 1);
      h_nnneighbour_temp[78 + b * 6 + 3].type._type = getType(h_nnneighbour_temp[78 + b * 6 + 3].id);

      h_nnneighbour_temp[78 + b * 6 + 4].id = lattice_list->getId(x4, y4 + 1, z4);
      h_nnneighbour_temp[78 + b * 6 + 4].type._type = getType(h_nnneighbour_temp[78 + b * 6 + 4].id);

      h_nnneighbour_temp[78 + b * 6 + 5].id = lattice_list->getId(x4 + 2, y4, z4);
      h_nnneighbour_temp[78 + b * 6 + 5].type._type = getType(h_nnneighbour_temp[78 + b * 6 + 5].id);

    }
    if(sect % 2 == 0) {
      if(h_event_id < vac_count_adv) {
        HANDLE_HIP(hipMemcpy(d_vacancy_adv + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        HANDLE_HIP(hipMemcpy(d_nnneighbour_adv + NN_TOTAL * h_event_id, h_nnneighbour_temp, NN_TOTAL * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
      }else {
        HANDLE_HIP(hipMemcpy(d_vacancy_aft + h_event_id - vac_count_adv, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        HANDLE_HIP(hipMemcpy(d_nnneighbour_aft + NN_TOTAL * (h_event_id - vac_count_adv), h_nnneighbour_temp, NN_TOTAL * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
      }
    }else {
      if(h_event_id < vac_count_buf_adv) {
        HANDLE_HIP(hipMemcpy(d_vacancy_buf_adv + h_event_id, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        HANDLE_HIP(hipMemcpy(d_nnneighbour_buf_adv + NN_TOTAL * h_event_id, h_nnneighbour_temp, NN_TOTAL * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
      }else {
        HANDLE_HIP(hipMemcpy(d_vacancy_buf_aft + h_event_id - vac_count_buf_adv, &h_vacancy_temp, sizeof(dev_Vacancy), hipMemcpyHostToDevice));
        HANDLE_HIP(hipMemcpy(d_nnneighbour_buf_aft + NN_TOTAL * (h_event_id - vac_count_buf_adv), h_nnneighbour_temp, NN_TOTAL * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
      }
    }
  }
    
  if(sect % 2 == 0) {

    blockSize = 64;
    numBlocks = (vac_count_adv + vac_count_aft + blockSize - 1) / blockSize;

    // 调用GPU内核函数
    // hipLaunchKernelGGL(updatennLatti, numBlocks, blockSize, 0, 0, 
    //                    d_event,
    //                    d_nnneighbour_adv,
    //                    d_nnneighbour_aft,
    //                    vac_count_adv,
    //                    vac_count_adv + vac_count_aft);
    // 算子融合，上面这个更新信息的 kernel 函数合并到计算 Kernel 里 
    hipLaunchKernelGGL(calcRatesGPU, numBlocks, blockSize, 0, 0, 
                       d_event,
                       d_vacancy_adv,
                       d_vacancy_aft,
                       d_nnneighbour_adv, 
                       d_nnneighbour_aft,
                       d_config, 
                       d_sum_rates, 
                       vac_count_adv,
                       vac_count_adv + vac_count_aft);
  }else {
    blockSize = 64;
    numBlocks = (vac_count_buf_adv + vac_count_buf_aft + blockSize - 1) / blockSize;

    // 调用GPU内核函数
    // hipLaunchKernelGGL(updatennLatti, numBlocks, blockSize, 0, 0, 
    //                    d_event,
    //                    d_nnneighbour_buf_adv,
    //                    d_nnneighbour_buf_aft,
    //                    vac_count_buf_adv,
    //                    vac_count_buf_adv + vac_count_buf_aft);
    // 算子融合，上面这个更新信息的 kernel 函数合并到计算 Kernel 里 
    hipLaunchKernelGGL(calcRatesGPU, numBlocks, blockSize, 0, 0,
                       d_event,
                       d_vacancy_buf_adv,
                       d_vacancy_buf_aft,
                       d_nnneighbour_buf_adv, 
                       d_nnneighbour_buf_aft,
                       d_config, 
                       d_sum_rates, 
                       vac_count_buf_adv,
                       vac_count_buf_adv + vac_count_buf_aft);
  }
  CHECK_FOR_ERROR();

  *sum_rates = 0;

  HANDLE_HIP(hipMemcpy(sum_rates, d_sum_rates, sizeof(double), hipMemcpyDeviceToHost));
  return *sum_rates;
}

void calculateNextAdvRegion(int sect){
  using namespace gpuSimulate;
  if(sect % 2 == 0){
    if(vac_count_buf_adv > 0) {
      blockSize = 64;
      numBlocks = (vac_count_buf_adv + blockSize - 1) / blockSize;
      // 调用GPU内核函数
      hipLaunchKernelGGL(calcRatesGPUAdv, numBlocks, blockSize, 0, stream_buffer[3], 
                         d_vacancy_buf_adv,
                         d_nnneighbour_buf_adv, 
                         d_config, 
                         d_sum_rates_adv, 
                         vac_count_buf_adv);
    }
  }else {
    if(vac_count_adv > 0) {
      blockSize = 64;
      numBlocks = (vac_count_adv + blockSize - 1) / blockSize;
      // 调用GPU内核函数
      hipLaunchKernelGGL(calcRatesGPUAdv, numBlocks, blockSize, 0, stream[3], 
                         d_vacancy_adv,
                         d_nnneighbour_adv, 
                         d_config, 
                         d_sum_rates_adv, 
                         vac_count_adv);
    }
  }
}

void sector_final(int sect){
  using namespace gpuSimulate;
  // kiwi::logs::v(" ", " vac_count_adv is : {} vac_count_aft is : {} vac_valid_count is : {} vac_count_buf_adv is : {} vac_count_buf_aft is : {} vac_valid_buf_count is : {}.\n", 
  //             vac_count_adv, vac_count_aft, vac_valid_count, vac_count_buf_adv, vac_count_buf_aft, vac_valid_buf_count);
  // kiwi::logs::v(" ", " sector_final !!!!!!!!!!! .\n");
  if(sect % 2 == 0) {
    vac_count_adv = 0;
    vac_count_aft = 0;
    vac_valid_count = 0;
  }else {
    vac_count_buf_adv = 0;
    vac_count_buf_aft = 0;
    vac_valid_buf_count = 0;
  }

}

void initialize_gpu(int process_rank){
    int gpu_count;
    HANDLE_HIP(hipGetDeviceCount(&gpu_count));
    HANDLE_HIP(hipSetDevice(process_rank % gpu_count));
    HANDLE_HIP(hipDeviceReset());
}

void gpu_prepare(const double v, const double T, LatticesList *p_list, comm::ColoredDomain *_p_domain){
  using namespace gpuSimulate;

  lattice_list = p_list;
  p_domain = _p_domain;
  double h_config[2];
  h_config[0] = v;   // v
  h_config[1] = T;   // T

  // HANDLE_HIP(hipMalloc((void**)&d_size, 2 * sizeof(long)));
  HANDLE_HIP(hipMalloc((void**)&d_config, 2 * sizeof(double)));

  //HANDLE_HIP(hipMemcpy(d_size, h_size, 2 * sizeof(long), hipMemcpyHostToDevice));
  HANDLE_HIP(hipMemcpy(d_config, h_config, 2 * sizeof(double), hipMemcpyHostToDevice));

  HANDLE_HIP(hipMalloc((void**)&d_sum_rates, sizeof(double)));
  HANDLE_HIP(hipMalloc((void**)&d_sum_rates_adv, sizeof(double)));
  HANDLE_HIP(hipMalloc((void**)&d_sum_rates_aft, sizeof(double)));
  HANDLE_HIP(hipMalloc((void**)&d_event, sizeof(dev_event)));
  HANDLE_HIP(hipMalloc((void**)&d_event_id, sizeof(_type_lattice_count)));

  hipHostMalloc((void**)&h_nnneighbour_temp, NN_TOTAL * sizeof(dev_nnLattice));
  
  vac_host_array_adv_size = lattice_list->vac_hash.size() / 2 + 1;
  vac_host_array_aft_size = lattice_list->vac_hash.size() / 2 + 1;
  vac_host_array_buf_adv_size = lattice_list->vac_hash.size() / 2 + 1;
  vac_host_array_buf_aft_size = lattice_list->vac_hash.size() / 2 + 1;

  HANDLE_HIP(hipMalloc((void**)&d_vacancy_adv, vac_host_array_adv_size * sizeof(dev_Vacancy)));
  HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_adv, vac_host_array_adv_size * NN_TOTAL * sizeof(dev_nnLattice)));

  HANDLE_HIP(hipMalloc((void**)&d_vacancy_aft, vac_host_array_aft_size * sizeof(dev_Vacancy)));
  HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_aft, vac_host_array_aft_size * NN_TOTAL * sizeof(dev_nnLattice)));

  HANDLE_HIP(hipMalloc((void**)&d_vacancy_buf_adv, vac_host_array_buf_adv_size * sizeof(dev_Vacancy)));
  HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_buf_adv, vac_host_array_buf_adv_size * NN_TOTAL * sizeof(dev_nnLattice)));

  HANDLE_HIP(hipMalloc((void**)&d_vacancy_buf_aft, vac_host_array_buf_aft_size * sizeof(dev_Vacancy)));
  HANDLE_HIP(hipMalloc((void**)&d_nnneighbour_buf_aft, vac_host_array_buf_aft_size * NN_TOTAL * sizeof(dev_nnLattice)));

  hipHostMalloc((void**)&sum_rates, sizeof(double));
  hipHostMalloc((void**)&sum_rates_adv, sizeof(double));
  hipHostMalloc((void**)&sum_rates_aft, sizeof(double));

  hipHostMalloc((void**)&h_vacancy_adv, vac_host_array_adv_size * sizeof(dev_Vacancy));
  hipHostMalloc((void**)&h_nnneighbour_adv, vac_host_array_adv_size * NN_TOTAL * sizeof(dev_nnLattice));

  hipHostMalloc((void**)&h_vacancy_aft, vac_host_array_aft_size * sizeof(dev_Vacancy));
  hipHostMalloc((void**)&h_nnneighbour_aft, vac_host_array_aft_size * NN_TOTAL * sizeof(dev_nnLattice));

  hipHostMalloc((void**)&h_vacancy_buf_adv, vac_host_array_buf_adv_size * sizeof(dev_Vacancy));
  hipHostMalloc((void**)&h_nnneighbour_buf_adv, vac_host_array_buf_adv_size * NN_TOTAL * sizeof(dev_nnLattice));

  hipHostMalloc((void**)&h_vacancy_buf_aft, vac_host_array_buf_aft_size * sizeof(dev_Vacancy));
  hipHostMalloc((void**)&h_nnneighbour_buf_aft, vac_host_array_buf_aft_size * NN_TOTAL * sizeof(dev_nnLattice));

  for(int i = 0; i < 4; i++){
    hipStreamCreate(&stream[i]);
    hipStreamCreate(&stream_buffer[i]);
  }
}

void gpu_final(){

    using namespace gpuSimulate;
    // 释放GPU内存
    HANDLE_HIP(hipFree(d_config));
    HANDLE_HIP(hipFree(d_sum_rates));
    HANDLE_HIP(hipFree(d_sum_rates_adv));
    HANDLE_HIP(hipFree(d_sum_rates_aft));
    HANDLE_HIP(hipFree(d_event));
    HANDLE_HIP(hipFree(d_event_id));

    HANDLE_HIP(hipFree(d_vacancy_adv));
    HANDLE_HIP(hipFree(d_nnneighbour_adv));
    HANDLE_HIP(hipFree(d_vacancy_aft));
    HANDLE_HIP(hipFree(d_nnneighbour_aft));
    HANDLE_HIP(hipFree(d_vacancy_buf_adv));
    HANDLE_HIP(hipFree(d_nnneighbour_buf_adv));
    HANDLE_HIP(hipFree(d_vacancy_buf_aft));
    HANDLE_HIP(hipFree(d_nnneighbour_buf_aft));

    hipHostFree(sum_rates);
    hipHostFree(sum_rates_adv);
    hipHostFree(sum_rates_aft);

    hipHostFree(h_nnneighbour_temp);

    hipHostFree(h_vacancy_adv);
    hipHostFree(h_nnneighbour_adv);
    hipHostFree(h_vacancy_aft);
    hipHostFree(h_nnneighbour_aft);
    hipHostFree(h_vacancy_buf_adv);
    hipHostFree(h_nnneighbour_buf_adv);
    hipHostFree(h_vacancy_buf_aft);
    hipHostFree(h_nnneighbour_buf_aft);

    for(int i = 0; i < 4; i++){
      hipStreamDestroy(stream[i]);
      hipStreamDestroy(stream_buffer[i]);
    }
}

void initInfo(_type_lattice_count vac_total, _type_lattice_id *vac_idArray,
              dev_Vacancy *h_vacancy, dev_nnLattice *h_nnneighbour) {

  using namespace gpuSimulate;
  // 准备可能用到的lattice
  
  for(_type_lattice_count i = 0; i < vac_total; i++){

    _type_lattice_id x = vac_idArray[i] % lattice_list->meta.size_x;
    _type_lattice_id y = (vac_idArray[i] / lattice_list->meta.size_x) % lattice_list->meta.size_y;
    _type_lattice_id z = vac_idArray[i] / (lattice_list->meta.size_x * lattice_list->meta.size_y);

    h_vacancy[i].id = vac_idArray[i];
    h_vacancy[i].valid = 1;

    if(x % 2 == 0){

      h_nnneighbour[i * NN_TOTAL].id = lattice_list->getId(x - 1, y - 1, z - 1);
      h_nnneighbour[i * NN_TOTAL].type._type = getType(h_nnneighbour[i * NN_TOTAL].id); // 126 = 8 + 6 + 64 + 48

      h_nnneighbour[i * NN_TOTAL + 1].id = lattice_list->getId(x - 1, y - 1, z);
      h_nnneighbour[i * NN_TOTAL + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 1].id);

      h_nnneighbour[i * NN_TOTAL + 2].id = lattice_list->getId(x - 1, y, z - 1);
      h_nnneighbour[i * NN_TOTAL + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 2].id);

      h_nnneighbour[i * NN_TOTAL + 3].id = lattice_list->getId(x - 1, y, z);
      h_nnneighbour[i * NN_TOTAL + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 3].id);

      h_nnneighbour[i * NN_TOTAL + 4].id = lattice_list->getId(x + 1, y - 1, z - 1);
      h_nnneighbour[i * NN_TOTAL + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 4].id);

      h_nnneighbour[i * NN_TOTAL + 5].id = lattice_list->getId(x + 1, y - 1, z);
      h_nnneighbour[i * NN_TOTAL + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 5].id);

      h_nnneighbour[i * NN_TOTAL + 6].id = lattice_list->getId(x + 1, y, z - 1);
      h_nnneighbour[i * NN_TOTAL + 6].type._type = getType(h_nnneighbour[i * NN_TOTAL + 6].id);

      h_nnneighbour[i * NN_TOTAL + 7].id = lattice_list->getId(x + 1, y, z);
      h_nnneighbour[i * NN_TOTAL + 7].type._type = getType(h_nnneighbour[i * NN_TOTAL + 7].id);

    }else{
        
      h_nnneighbour[i * NN_TOTAL].id = lattice_list->getId(x - 1, y, z);
      h_nnneighbour[i * NN_TOTAL].type._type = getType(h_nnneighbour[i * NN_TOTAL].id);

      h_nnneighbour[i * NN_TOTAL + 1].id = lattice_list->getId(x - 1, y, z + 1);
      h_nnneighbour[i * NN_TOTAL + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 1].id);

      h_nnneighbour[i * NN_TOTAL + 2].id = lattice_list->getId(x - 1, y + 1, z);
      h_nnneighbour[i * NN_TOTAL + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 2].id);

      h_nnneighbour[i * NN_TOTAL + 3].id = lattice_list->getId(x - 1, y + 1, z + 1);
      h_nnneighbour[i * NN_TOTAL + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 3].id);

      h_nnneighbour[i * NN_TOTAL + 4].id = lattice_list->getId(x + 1, y, z);
      h_nnneighbour[i * NN_TOTAL + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 4].id);

      h_nnneighbour[i * NN_TOTAL + 5].id = lattice_list->getId(x + 1, y, z + 1);
      h_nnneighbour[i * NN_TOTAL + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 5].id);

      h_nnneighbour[i * NN_TOTAL + 6].id = lattice_list->getId(x + 1, y + 1, z);
      h_nnneighbour[i * NN_TOTAL + 6].type._type = getType(h_nnneighbour[i * NN_TOTAL + 6].id);

      h_nnneighbour[i * NN_TOTAL + 7].id = lattice_list->getId(x + 1, y + 1, z + 1);
      h_nnneighbour[i * NN_TOTAL + 7].type._type = getType(h_nnneighbour[i * NN_TOTAL + 7].id);

    }

    h_nnneighbour[i * NN_TOTAL + 8].id = lattice_list->getId(x - 2, y, z);
    h_nnneighbour[i * NN_TOTAL + 8].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8].id);

    h_nnneighbour[i * NN_TOTAL + 8 + 1].id = lattice_list->getId(x, y - 1, z);
    h_nnneighbour[i * NN_TOTAL + 8 + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8 + 1].id);

    h_nnneighbour[i * NN_TOTAL + 8 + 2].id = lattice_list->getId(x, y, z - 1);
    h_nnneighbour[i * NN_TOTAL + 8 + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8 + 2].id);

    h_nnneighbour[i * NN_TOTAL + 8 + 3].id = lattice_list->getId(x, y, z + 1);
    h_nnneighbour[i * NN_TOTAL + 8 + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8 + 3].id);

    h_nnneighbour[i * NN_TOTAL + 8 + 4].id = lattice_list->getId(x, y + 1, z);
    h_nnneighbour[i * NN_TOTAL + 8 + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8 + 4].id);

    h_nnneighbour[i * NN_TOTAL + 8 + 5].id = lattice_list->getId(x + 2, y, z);
    h_nnneighbour[i * NN_TOTAL + 8 + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 8 + 5].id);

    for(int b = 0; b < 8; b++){
      _type_lattice_id x2 = h_nnneighbour[i * NN_TOTAL + b].id % lattice_list->meta.size_x;
      _type_lattice_id y2 = (h_nnneighbour[i * NN_TOTAL + b].id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
      _type_lattice_id z2 = h_nnneighbour[i * NN_TOTAL + b].id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

      if(x2 % 2 == 0){
          
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8].id = lattice_list->getId(x2 - 1, y2 - 1, z2 - 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8].id); // 14 = 6 + 8

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].id = lattice_list->getId(x2 - 1, y2 - 1, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].id = lattice_list->getId(x2 - 1, y2, z2 - 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].id = lattice_list->getId(x2 - 1, y2, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].id = lattice_list->getId(x2 + 1, y2 - 1, z2 - 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].id = lattice_list->getId(x2 + 1, y2 - 1, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].id = lattice_list->getId(x2 + 1, y2, z2 - 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].id = lattice_list->getId(x2 + 1, y2, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].id);

      }else{

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8].id = lattice_list->getId(x2 - 1, y2, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].id = lattice_list->getId(x2 - 1, y2, z2 + 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 1].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].id = lattice_list->getId(x2 - 1, y2 + 1, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 2].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].id = lattice_list->getId(x2 - 1, y2 + 1, z2 + 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 3].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].id = lattice_list->getId(x2 + 1, y2, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 4].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].id = lattice_list->getId(x2 + 1, y2, z2 + 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 5].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].id = lattice_list->getId(x2 + 1, y2 + 1, z2);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 6].id);

        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].id = lattice_list->getId(x2 + 1, y2 + 1, z2 + 1);
        h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].type._type = getType(h_nnneighbour[i * NN_TOTAL + 14 + b * 8 + 7].id);

      }

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6].id = lattice_list->getId(x2 - 2, y2, z2);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6].id); // 78 = 64 + 14

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 1].id = lattice_list->getId(x2, y2 - 1, z2);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 1].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 1].id);

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 2].id = lattice_list->getId(x2, y2, z2 - 1);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 2].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 2].id);

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 3].id = lattice_list->getId(x2, y2, z2 + 1);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 3].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 3].id);

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 4].id = lattice_list->getId(x2, y2 + 1, z2);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 4].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 4].id);

      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 5].id = lattice_list->getId(x2 + 2, y2, z2);
      h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 5].type._type = getType(h_nnneighbour[i * NN_TOTAL + 78 + b * 6 + 5].id);

    }
  }
}

LatticeTypes::lat_type getType(_type_lattice_id latti_id){
  using namespace gpuSimulate;

  auto it_vac = lattice_list->vac_hash.find(latti_id);
  auto it_cu = lattice_list->cu_hash.find(latti_id);
  auto it_mn = lattice_list->mn_hash.find(latti_id);
  auto it_ni = lattice_list->ni_hash.find(latti_id);
  if (it_vac != lattice_list->vac_hash.end()) {
    return LatticeTypes::V;
  }else if (it_cu != lattice_list->cu_hash.end()){
    return LatticeTypes::Cu;
  }else if (it_mn != lattice_list->mn_hash.end()){
    return LatticeTypes::Mn;
  }else if (it_ni != lattice_list->ni_hash.end()){
    return LatticeTypes::Ni;
  }else{
    return LatticeTypes::Fe;
  }
}

void exchange_lattices(dev_event h_event){
  using namespace gpuSimulate;

  auto it_vac = lattice_list->vac_hash.find(h_event.from_id);
  if (it_vac != lattice_list->vac_hash.end()) {
    lattice_list->vac_hash.erase(it_vac);
    lattice_list->vac_hash.emplace(std::make_pair(h_event.to_id, VacancyHash{}));
  }else{
    assert(false);
  }

  switch (h_event.to_type._type)
  {
  case LatticeTypes::Cu:
    lattice_list->cu_hash.erase(h_event.to_id);
    lattice_list->cu_hash.emplace(h_event.from_id);
    break;
  case LatticeTypes::Mn:
    lattice_list->mn_hash.erase(h_event.to_id);
    lattice_list->mn_hash.emplace(h_event.from_id);
    break;
  case LatticeTypes::Ni:
    lattice_list->ni_hash.erase(h_event.to_id);
    lattice_list->ni_hash.emplace(h_event.from_id);
    break;
  default:
    //
    break;
  }
}

void add_ghost(std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
               _type_lattice_id x, _type_lattice_id y, _type_lattice_id z, Lattice lat_to, const unsigned int sector_id){

  using namespace gpuSimulate;
      
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

void add_surface(_type_lattice_id surface_id,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z) {
  using namespace gpuSimulate;
  _type_lattice_id surface_x = surface_id % lattice_list->meta.size_x;
  _type_lattice_id surface_y = (surface_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  _type_lattice_id surface_z = surface_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  // todo the regions can be static.
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  // 应该可以不循环 8 次的 ？ 如果发送区没有重叠的话
  for (int sect = 0; sect < 8; sect++) {
    for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
      send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
    }
    for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
      for (auto &r : send_regions[d]) {
        _type_lattice_size x_low = 2 * r.x_low;
        _type_lattice_size y_low = r.y_low;
        _type_lattice_size z_low = r.z_low;
        _type_lattice_size x_high = 2 * r.x_high;
        _type_lattice_size y_high = r.y_high;
        _type_lattice_size z_high = r.z_high;
    
        if(x_low <= surface_x && surface_x < x_high && y_low <= surface_y && surface_y < y_high && z_low <= surface_z && surface_z < z_high){
          if(d == 0){
            auto it_from = exchange_surface_x[sect].find(surface_id);
            if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(surface_id);
          }else if(d == 1) {
            auto it_from = exchange_surface_y[sect].find(surface_id);
            if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(surface_id);
          }else {
            auto it_from = exchange_surface_z[sect].find(surface_id);
            if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(surface_id);
          }
        }
      }
    }
  }
}