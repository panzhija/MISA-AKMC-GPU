//
// Created by genshen on 2018-12-12.
//

#include <cassert>

#include "env.h"
#include "kmc.h"
#include "rate/itl_rates_solver.h"
#include "rate/vacancy_rates_solver.h"
#include "recombine.h"
#include "utils/random/random.h"
#include <iostream>
#include <logs/logs.h>
#include "../gpu/gpu_simulate.h"
#include "utils/simulation_domain.h"
#include <map>
#include <vector>
#include <algorithm>
ABVIModel::ABVIModel(Box *box, double v, double T) : box(box), v(v), T(T) {}

void ABVIModel::transFirstsectLocalToGpu(const comm::Region<comm::_type_lattice_size> region) {

  _type_lattice_id *vac_idArray_adv;
  vac_idArray_adv = (_type_lattice_id *)malloc(box->lattice_list->vac_hash.size() * sizeof(_type_lattice_id));
  _type_lattice_count vac_adv = 0;

  _type_lattice_size x_adv_low = region.x_low + p_domain->lattice_size_ghost[0] + 3;
  _type_lattice_size y_adv_low = region.y_low + p_domain->lattice_size_ghost[1] + 3;
  _type_lattice_size z_adv_low = region.z_low + p_domain->lattice_size_ghost[2] + 3;
  _type_lattice_size x_adv_high = region.x_high - p_domain->lattice_size_ghost[0] - 3;
  _type_lattice_size y_adv_high = region.y_high - p_domain->lattice_size_ghost[1] - 3;
  _type_lattice_size z_adv_high = region.z_high - p_domain->lattice_size_ghost[2] - 3;

  assert(x_adv_low < x_adv_high);
  assert(y_adv_low < y_adv_high);
  assert(z_adv_low < z_adv_high);
  for (const auto& pair : box->lattice_list->vac_hash) {
    _type_lattice_size x = pair.first % box->lattice_list->meta.size_x;
    _type_lattice_size y = (pair.first / box->lattice_list->meta.size_x) % box->lattice_list->meta.size_y;
    _type_lattice_size z = pair.first / (box->lattice_list->meta.size_x * box->lattice_list->meta.size_y);

    if(2 * x_adv_low <= x && x < x_adv_high && y_adv_low <= y && y < y_adv_high && z_adv_low <= z && z < z_adv_high){
      vac_idArray_adv[vac_adv++] = pair.first;
      // kiwi::logs::v(" ", " id is : {} x is : {} y is : {} z is : {}.\n", pair.first, x, y, z);
    }
  }
  // kiwi::logs::v(" ", " first vac_adv count is : {}.\n", vac_adv);
  firstSectLocalToGpu(vac_idArray_adv, vac_adv);
}

// 这个函数的主要作用是计算指定仿真区域内的某种类型的迁移速率（transition rates）。
// 在材料模拟或分子动力学中，迁移速率通常用于描述原子、分子或其他粒子之间的运动、迁移或反应速度。具体来说，这个函数执行以下操作：
// 初始化一些变量和对象，如ItlRatesSolver和VacRatesSolver，这些对象可能用于计算迁移速率。
// 使用嵌套循环遍历指定仿真区域内的所有格点。这些嵌套循环涵盖三个维度：x、y 和 z。
// 对于每个格点，检查格点上的晶格类型（Lattice type），可能是"dumbbell"或"vacancy"，分别代表"双原子晶格"和"空位"。根据不同的晶格类型，选择不同的处理路径。
// 对于"dumbbell"类型的晶格，执行一系列操作来计算与相邻格点之间的迁移速率。这包括获取与当前晶格相邻的格点信息、更新迁移速率等。
// 迁移速率计算似乎涉及到一个lambda函数，该函数返回从当前晶格到相邻晶格的迁移速率，并将其累加到sum_rates中。
// 对于"vacancy"类型的晶格，执行类似的操作，计算空位的迁移速率，并将其累加到sum_rates中。
// 最后，通过调用defectGenRate()函数，将生成的缺陷（defect）的生成速率添加到sum_rates中。
// 返回sum_rates，即在整个仿真区域内计算得到的总迁移速率。
// 总的来说，用于计算指定仿真区域内不同晶格类型的迁移速率，并将它们累加到一起以获得总迁移速率。

// _type_rate ABVIModel::calcRates_GPU(const comm::Region<comm::_type_lattice_size> region) {
//   return calcRatesGPU(region, box, v, T);
// }

_type_rate ABVIModel::calcRates(const comm::Region<comm::_type_lattice_size> region) {
  _type_rate sum_rates = 0;
  VacRatesSolver vac_rate(*(box->lattice_list), v, T);

  for (_type_lattice_id z = region.z_low; z < region.z_high; z++) {
    for (_type_lattice_id y = region.y_low; y < region.y_high; y++) {
      for (_type_lattice_id x = 2 * region.x_low; x < 2 * region.x_high; x++) {
        _type_lattice_id latti_id = box->lattice_list->getId(x, y, z);
        auto it_vac = box->lattice_list->vac_hash.find(latti_id);
        if (it_vac != box->lattice_list->vac_hash.end()) { // vacancy
#ifdef KMC_DEBUG_MODE
          // assert(box->va_list->mp.find(lattice.getId()) != box->va_list->mp.end());
#endif
          // _type_lattice_size x = latti_id % box->lattice_list->meta.size_x;
          // _type_lattice_size y = (latti_id / box->lattice_list->meta.size_x) % box->lattice_list->meta.size_y;
          // _type_lattice_size z = latti_id / (box->lattice_list->meta.size_x * box->lattice_list->meta.size_y);
          // kiwi::logs::v(" ", " id is : {} x is : {} y is : {} z is : {}.\n", latti_id, x, y, z);
          VacancyHash &vacancyhash = box->lattice_list->vac_hash.at(latti_id);
          Lattice lat_list[LatticesList::MAX_1NN]; // todo new array many times.
          box->lattice_list->get1nn2(x, y, z, lat_list);
          _type_neighbour_status nei_status = box->lattice_list->get1nnStatus(x, y, z);
          //kiwi::logs::v(" ", " nei_status is : {} .\n", nei_status);
          vacancyhash.beforeRatesUpdate2(lat_list, nei_status);
          // for(int i = 0; i < 8; i++){
          //   kiwi::logs::v(" ", " {} is : {} .\n", i, lat_list[i].type._type);
          // }
          Lattice source_latti;
          source_latti.id = latti_id;
          source_latti.type = LatticeTypes{LatticeTypes::V};
          //kiwi::logs::v(" ", " avail_trans_dir is : {} .\n", vacancy.avail_trans_dir);
          vacancyhash.updateRates2(source_latti, lat_list, nei_status,
                              [&source_latti, &vac_rate, &sum_rates](Lattice *lat_nei,
                                                                const LatticeTypes::lat_type ghost_atom,
                                                                const _type_dir_id _1nn_offset) -> _type_rate {
                                _type_rate rate = vac_rate.rate(source_latti, *lat_nei, ghost_atom, _1nn_offset);
                                sum_rates += rate; // add this rate to sum
                                //kiwi::logs::v(" ", " rates is : {} .\n", rate);
                                return rate;
                              });
          
        }
        // there is no transition rate for single atom.
      }
    }
  }
  sum_rates += defectGenRate();
  return sum_rates;
}


_type_rate ABVIModel::calcRatesGPU(const comm::Region<comm::_type_lattice_size> region, bool *sector_first, const comm::Region<comm::_type_lattice_size> region_next, int sect, int is_last_step) {
  // for(std::map<_type_lattice_id, Vacancy>::iterator iter = box->va_list->mp.begin(); iter != box->va_list->mp.end(); ++iter){
  //   std::cout << "rank is" << SimulationDomain::comm_sim_pro.own_rank << "vac_map->first: " << iter->first << std::endl;
  // }
  if(*sector_first){
    _type_lattice_id *vac_idArray_aft;
    vac_idArray_aft = (_type_lattice_id *)malloc(box->lattice_list->vac_hash.size() * sizeof(_type_lattice_id));
    _type_lattice_count vac_aft = 0;

    _type_lattice_id *vac_idArray_buf_adv;
    vac_idArray_buf_adv = (_type_lattice_id *)malloc(box->lattice_list->vac_hash.size() * sizeof(_type_lattice_id));
    _type_lattice_count vac_adv = 0;

    _type_lattice_size x_adv_low = region.x_low + p_domain->lattice_size_ghost[0] + 3;
    _type_lattice_size y_adv_low = region.y_low + p_domain->lattice_size_ghost[1] + 3;
    _type_lattice_size z_adv_low = region.z_low + p_domain->lattice_size_ghost[2] + 3;
    _type_lattice_size x_adv_high = region.x_high - p_domain->lattice_size_ghost[0] - 3;
    _type_lattice_size y_adv_high = region.y_high - p_domain->lattice_size_ghost[1] - 3;
    _type_lattice_size z_adv_high = region.z_high - p_domain->lattice_size_ghost[2] - 3;
    assert(x_adv_low < x_adv_high);
    assert(y_adv_low < y_adv_high);
    assert(z_adv_low < z_adv_high);
    
    _type_lattice_size x_adv_buf_low = region_next.x_low + p_domain->lattice_size_ghost[0] + 3;
    _type_lattice_size y_adv_buf_low = region_next.y_low + p_domain->lattice_size_ghost[1] + 3;
    _type_lattice_size z_adv_buf_low = region_next.z_low + p_domain->lattice_size_ghost[2] + 3;
    _type_lattice_size x_adv_buf_high = region_next.x_high - p_domain->lattice_size_ghost[0] - 3;
    _type_lattice_size y_adv_buf_high = region_next.y_high - p_domain->lattice_size_ghost[1] - 3;
    _type_lattice_size z_adv_buf_high = region_next.z_high - p_domain->lattice_size_ghost[2] - 3;
    assert(x_adv_buf_low < x_adv_buf_high);
    assert(y_adv_buf_low < y_adv_buf_high);
    assert(z_adv_buf_low < z_adv_buf_high);

    bool is_adv;
  
    for (const auto& pair : box->lattice_list->vac_hash) {
      _type_lattice_size x = pair.first % box->lattice_list->meta.size_x;
      _type_lattice_size y = (pair.first / box->lattice_list->meta.size_x) % box->lattice_list->meta.size_y;
      _type_lattice_size z = pair.first / (box->lattice_list->meta.size_x * box->lattice_list->meta.size_y);
      is_adv = 2 * x_adv_low <= x && x < x_adv_high && y_adv_low <= y && y < y_adv_high && z_adv_low <= z && z < z_adv_high;
      // 传输当前扇区剩下的空位
      if(2 * region.x_low <= x && x < 2 * region.x_high && region.y_low <= y && y < region.y_high && region.z_low <= z && z < region.z_high && !(is_adv)){
        vac_idArray_aft[vac_aft++] = pair.first;
        // kiwi::logs::v(" ", " id is : {} x is : {} y is : {} z is : {}.\n", pair.first, x, y, z);
      }
      // 若当前的扇区不是最后一次迭代的最后一个扇区，则预传下一个扇区的部分空位
      if(is_last_step != 1 || sect != 7) {
        if(2 * x_adv_buf_low <= x && x < x_adv_buf_high && y_adv_buf_low <= y && y < y_adv_buf_high && z_adv_buf_low <= z && z < z_adv_buf_high) {
          vac_idArray_buf_adv[vac_adv++] = pair.first;
          // kiwi::logs::v(" ", " id is : {} x is : {} y is : {} z is : {}.\n", pair.first, x, y, z);
        }
      }
    }
    *sector_first = false;
    // std::sort(vac_idArray, vac_idArray + vac_total);
    // kiwi::logs::v(" ", " vac_aft count is : {}.\n", vac_aft);
    // kiwi::logs::v(" ", " next_vac_adv count is : {}.\n", vac_adv);
    return first_calculate_GPU(vac_idArray_aft, vac_aft, vac_idArray_buf_adv, vac_adv, sect);
  }else{
    return calculate_GPU(region, sect);
  }
}

// _type_rate ABVIModel::calcRatesGPU(const comm::Region<comm::_type_lattice_size> region, bool *sector_first) {
//   // for(std::map<_type_lattice_id, Vacancy>::iterator iter = box->va_list->mp.begin(); iter != box->va_list->mp.end(); ++iter){
//   //   std::cout << "rank is" << SimulationDomain::comm_sim_pro.own_rank << "vac_map->first: " << iter->first << std::endl;
//   // }
//     _type_lattice_id vac_idArray[box->va_list->mp.size()];
//     int vac_total = 0;
//     for (int z = region.z_low; z < region.z_high; z++) {
//       for (int y = region.y_low; y < region.y_high; y++) {
//         for (int x = 2 * region.x_low; x < 2 * region.x_high; x++) {
//           Lattice &lattice = box->lattice_list->getLat(x, y, z); 
//           if (lattice.type.isVacancy()) {
//             vac_idArray[vac_total] = lattice.getId();
//             vac_total ++;
//           }
//         }
//       }
//     }
//     // kiwi::logs::v(" ", "vac_total is : {} .\n", vac_total);
//     if(vac_total == 0){
//       return 0;
//     }else{
//       //std::cout << "vac_total is : " << vac_total<< std::endl;
//       return first_calculate_GPU2(box->lattice_list, vac_idArray, vac_total, box);
//     }
// }

_type_rate ABVIModel::defectGenRate() { return env::global_env.defect_gen_rate; }

void ABVIModel::selectAndPerformOnGPU(const _type_rate rate, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) {
  selectAndPerformEventGPU(rate, rank, step, sect, exchange_ghost, sector_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
}

// lat_region region：表示某个区域的参数。
// _type_rate excepted_rate：期望的速率。取值在[0, sum_rates)
// _type_rate sum_rates：总速率。
// 函数内部操作：
// rate_accumulator：初始化一个机率累加器，用于追踪已累加的机率。
// selected_event：创建一个 event::SelectedEvent 对象，用于存储所选事件的信息，初始事件类型为 event::DefectGen（缺陷生成）。
// 通过三个嵌套的循环遍历指定区域内的晶格：
// 外循环遍历 z 轴范围。
// 中循环遍历 y 轴范围。
// 内循环遍历 x 轴范围。
// 对于每个晶格，检查晶格的类型：
// 如果是 dumbbell 类型的晶格，进入相应的分支：
// 遍历 itl_ref 对象的机率数组，累加机率到 rate_accumulator 中。
// 如果 rate_accumulator 大于 excepted_rate，则表示找到了事件：
// 获取与 rate_index 相关的相邻晶格 _1nn_list。
// 计算事件类型为 event::DumbbellTrans，并设置相关的事件信息，包括起始晶格、目标晶格、目标标签和旋转方向。
// 使用 goto 标签 EVENT_FOUND 跳出循环。
// 如果是 vacancy 类型的晶格，进入相应的分支：
// 遍历 vacancy 对象的速率数组，累加速率到 rate_accumulator 中。
// 如果 rate_accumulator 大于 excepted_rate，则表示找到了事件：
// 获取与 rate_index 相关的相邻晶格 _1nn_list。
// 计算事件类型为 event::VacancyTrans，并设置相关的事件信息，包括起始晶格、目标晶格、目标标签。
// 使用 goto 标签 EVENT_FOUND 跳出循环。
// 如果以上两种类型都不是，表示事件未找到。
// 最后，通过 EVENT_FOUND 跳转标签，返回所选的事件。
// 该函数的目的是在指定的区域内选择一个事件，该事件的速率满足 excepted_rate，并返回一个包含事件信息的 event::SelectedEvent 对象。
// 这个函数通过遍历晶格、累加速率、检查类型等方式来选择事件，并在找到事件后使用 goto 跳转标签返回所选的事件。
// 如果没有找到满足条件的事件，函数将返回一个默认的 event::DefectGen 事件。

event::SelectedEvent ABVIModel::select(const lat_region region, const _type_rate excepted_rate,
                                       const _type_rate sum_rates) {
  _type_rate rate_accumulator = 0.0;
  event::SelectedEvent selected_event{event::DefectGen, 0, 0}; // default event is defect generation.
  for (_type_lattice_size z = region.z_low; z < region.z_high; z++) {
    for (_type_lattice_size y = region.y_low; y < region.y_high; y++) {
      for (_type_lattice_size x = 2 * region.x_low; x < 2 * region.x_high; x++) {
        // Lattice &lattice = box->lattice_list->getLat(x, y, z);
        // filter defect type
//         if (lattice.type.isDumbbell()) { // dumbbell
//           const Itl &itl_ref = box->itl_list->mp.at(lattice.getId());
//           for (int rate_index = 0; rate_index < Itl::RATES_SIZE; rate_index++) {
//             rate_accumulator += itl_ref.rates[rate_index];
//             if (rate_accumulator > excepted_rate) {
//               Lattice *_1nn_list[LatticesList::MAX_1NN] = {nullptr};
//               box->lattice_list->get1nn(lattice.getId(), _1nn_list);
//               // convert rate index to 1nn list array index.
//               const orientation ori = box->itl_list->mp[lattice.getId()].orient;
//               const _type_dir_id _1nn_tag = Itl::get1nnIdByRatesIndex(rate_index, ori.availTransDirs());
// #ifdef KMC_DEBUG_MODE
//               assert(_1nn_list[_1nn_tag] != nullptr);
// #endif
//               selected_event.event_type = event::DumbbellTrans;
//               selected_event.from_id = lattice.getId();
//               selected_event.to_id = box->lattice_list->meta.getIdBy1nnOffset(lattice.getId(), _1nn_tag);
//               selected_event.target_tag = _1nn_tag;
//               selected_event.rotate_direction = rate_index % 2 != 0;
//               goto EVENT_FOUND; // event found
//             }
//           }
        _type_lattice_size latti_id = box->lattice_list->getId(x, y, z);
        auto it_vac = box->lattice_list->vac_hash.find(latti_id);
        if (it_vac != box->lattice_list->vac_hash.end()) { // vacancy
          const VacancyHash &vacancy = box->lattice_list->vac_hash.at(latti_id);
          for (int rate_index = 0; rate_index < Vacancy::RATES_SIZE; rate_index++) {
            rate_accumulator += vacancy.nn_rates[rate_index];
            //kiwi::logs::v(" ", " {} rates is : {} .\n", rate_index, vacancy.nn_rates[rate_index]);
            if (rate_accumulator > excepted_rate) {
              // Lattice *_1nn_list[LatticesList::MAX_1NN] = {nullptr};
              // box->lattice_list->get1nn(lattice.getId(), _1nn_list);
#ifdef KMC_DEBUG_MODE
              assert(_1nn_list[rate_index] != nullptr);
#endif
              selected_event.event_type = event::VacancyTrans;
              selected_event.from_id = latti_id;
              selected_event.to_id = box->lattice_list->meta.getIdBy1nnOffset(latti_id, rate_index);
              // assert(selected_event.to_id == _1nn_list[rate_index]->getId());
              selected_event.target_tag = static_cast<_type_dir_id>(rate_index);
#ifdef KMC_DEBUG_MODE
              auto from_lat = box->lattice_list->getLat(selected_event.from_id);
              auto to_lat = box->lattice_list->getLat(selected_event.to_id);
              assert(from_lat.type.isVacancy());
              assert(to_lat.type.isAtom());
#endif
              goto EVENT_FOUND; // event found
            }
          }
        }
        // event not found
      }
    }
  }
EVENT_FOUND:
#ifdef KMC_DEBUG_MODE
  //    todo assert total rates == rate_accumulator + defect_gen_rate
#endif
  return selected_event;
}
// 使用 switch 语句根据所选事件的类型进行不同的操作。
// 对于不同类型的事件，分别执行以下操作：
// event::VacancyTrans 事件：
// 获取起始晶格 lat_from 和目标晶格 lat_to。
// 交换晶格类型，即将 lat_from 的类型与 lat_to 的类型互换。
// 更新空位列表 box->va_list，将 lat_from 移除并添加 lat_to。
// 如果存在事件监听器 p_event_listener，则调用 onVacancyTrans 回调通知事件发生。
// 进行重新组合（recombination）操作，创建 rec::RecList 对象，选择最小的重新组合事件并执行。
// event::DumbbellTrans 事件：
// 获取起始晶格 lat_from 和目标晶格 lat_to。
// 复制起始晶格的相关信息，并保存在 itl_copy 中。
// 通过计算得到跳跃原子 jump_atom，并更新晶格类型以实现原子的交换。
// 更新晶格的方向信息。
// 如果存在事件监听器 p_event_listener，则调用 onDumbbellTrans 回调通知事件发生。
// 进行重新组合（recombination）操作，创建 rec::RecList 对象，选择最小的重新组合事件并执行。
// event::DefectGen 事件：
// 通过随机选择一个起始晶格，搜索满足条件的晶格对 lat_1 和 lat_2。
// 随机选择一个晶格作为空位（V），另一个晶格作为 dumbbell，然后交换晶格类型。
// 如果存在事件监听器 p_event_listener，则调用 onDefectGenerate 回调通知事件发生。
// 更新晶格的方向信息和相应的列表。
// 该函数的目的是执行所选的事件，根据事件类型执行不同的操作，包括交换晶格类型、更新空位列表、通知事件监听器以及进行重新组合操作等。
// 这个函数的实现根据不同事件类型的特点，进行了相应的处理。

// event::VacancyTrans 事件：
// 这个事件表示一个空位（Vacancy）向一个晶格（Atom）的转移。
// 首先，获取起始晶格 lat_from 和目标晶格 lat_to。
// 接下来，交换这两个晶格的类型，即将 lat_from 的类型和 lat_to 的类型互换。
// 更新空位列表 box->va_list，将 lat_from 移除，并将 lat_to 添加到空位列表中，表示空位移动到了新的位置。
// 如果存在事件监听器 p_event_listener，则调用 onVacancyTrans 回调通知事件发生，传递有关事件的相关信息。
// 最后，执行重新组合（recombination）操作，创建 rec::RecList 对象，选择最小的重新组合事件并执行。这个步骤可能导致一些额外的晶格类型的变化。
// event::DumbbellTrans 事件：
// 这个事件表示一个 dumbbell（双原子晶格）向一个晶格（Atom）的转移。
// 类似于上述步骤，首先获取起始晶格 lat_from 和目标晶格 lat_to。
// 然后，复制起始晶格的相关信息，并保存在 itl_copy 中，以便稍后使用。
// 通过计算得到跳跃原子 jump_atom，并更新晶格类型以实现原子的交换。例如，可以将 lat_from 中的 jump_atom 移动到 lat_to 中，同时将 lat_from 中的其他原子移动到 lat_to 中。
// 更新晶格的方向信息，以反映跳跃动作。
// 如果存在事件监听器 p_event_listener，则调用 onDumbbellTrans 回调通知事件发生，传递与事件相关的信息，如起始晶格类型、目标晶格类型等。
// 最后，执行重新组合（recombination）操作，创建 rec::RecList 对象，选择最小的重新组合事件并执行。这个步骤可能导致一些额外的晶格类型的变化。
// event::DefectGen 事件：
// 这个事件表示在晶格中生成一个缺陷。
// 通过随机选择一个起始晶格 lat_1，然后搜索附近的晶格以找到满足条件的 lat_2。
// 随机选择一个晶格作为空位（V），另一个晶格作为 dumbbell，然后交换它们的类型，以模拟缺陷的生成。
// 如果存在事件监听器 p_event_listener，则调用 onDefectGenerate 回调通知事件发生，传递与事件相关的信息，如起始晶格类型、目标晶格类型等。
// 更新晶格的方向信息，通常使用随机生成的方向。
// 最后，将生成的晶格和缺陷添加到相应的列表中，以反映新的晶格状态。
// 总之，perform 函数的目的是在不同的事件类型下执行相应的操作，以模拟晶格中事件的发生和影响。
// 这些事件可能涉及晶格类型的交换、空位的生成或移动、重新组合等，这些操作是基于模型和物理规律来执行的，以模拟晶格的动态行为。

void ABVIModel::perform(const event::SelectedEvent selected, const lat_region region, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) {
  // 获取转移前lattice的信息 引用！！
  // kiwi::logs::v(" ", " from_id is : {} to_id is : {} .\n", selected.from_id, selected.to_id);
  if(selected.event_type == event::VacancyTrans){
    auto it_vac = box->lattice_list->vac_hash.find(selected.from_id);
    if (it_vac != box->lattice_list->vac_hash.end()) {
      box->lattice_list->vac_hash.erase(it_vac);
      box->lattice_list->vac_hash.emplace(std::make_pair(selected.to_id, VacancyHash{}));
    }else{
      assert(false);
    }

    auto it_cu = box->lattice_list->cu_hash.find(selected.to_id);
    auto it_mn = box->lattice_list->mn_hash.find(selected.to_id);
    auto it_ni = box->lattice_list->ni_hash.find(selected.to_id);
    if (it_cu != box->lattice_list->cu_hash.end()){
      box->lattice_list->cu_hash.erase(it_cu);
      box->lattice_list->cu_hash.emplace(selected.from_id);
    }else if (it_mn != box->lattice_list->mn_hash.end()){
      box->lattice_list->mn_hash.erase(it_mn);
      box->lattice_list->mn_hash.emplace(selected.from_id);
    }else if (it_ni != box->lattice_list->ni_hash.end()){
      box->lattice_list->ni_hash.erase(it_ni);
      box->lattice_list->ni_hash.emplace(selected.from_id);
    }else{
      // Fe 不做操作
    }
    if (box->lattice_list->meta.isGhostLat(selected.to_id)){
      Lattice lat_to;
      lat_to.id = selected.to_id;
      lat_to.type = LatticeTypes{LatticeTypes::V};
      addExchange_ghost(lat_to, sector_id);
    }
    if (box->lattice_list->meta.isSurfaceLat(selected.from_id)) addExchange_surface(selected.from_id);
    if (box->lattice_list->meta.isSurfaceLat(selected.to_id)) addExchange_surface(selected.to_id);
  }else{
    assert(false);
  }
}

void ABVIModel::reindex(const lat_region region) {
  box->va_list->reindex(box->lattice_list, region);
  // todo reindex interval list
  // todo refresh counter
}
void ABVIModel::setEventListener(EventListener *p_listener) { p_event_listener = p_listener; }

void ABVIModel::setColoredDomain(comm::ColoredDomain *_p_domain) { p_domain = _p_domain; }

void ABVIModel::addExchange_ghost(Lattice lat_to, const unsigned int sector_id){
  _type_lattice_id x = lat_to.id % box->lattice_list->meta.size_x;
  _type_lattice_id y = (lat_to.id / box->lattice_list->meta.size_x) % box->lattice_list->meta.size_y;
  _type_lattice_id z = lat_to.id / (box->lattice_list->meta.size_x * box->lattice_list->meta.size_y);
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

void ABVIModel::addExchange_surface(_type_lattice_id surface_id) {
  _type_lattice_id surface_x = surface_id % box->lattice_list->meta.size_x;
  _type_lattice_id surface_y = (surface_id / box->lattice_list->meta.size_x) % box->lattice_list->meta.size_y;
  _type_lattice_id surface_z = surface_id / (box->lattice_list->meta.size_x * box->lattice_list->meta.size_y);
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  // todo the regions can be static.
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  _type_lattice_size x_low;
  _type_lattice_size y_low;
  _type_lattice_size z_low;
  _type_lattice_size x_high;
  _type_lattice_size y_high;
  _type_lattice_size z_high;
  for (int sect = 0; sect < 8; sect++) {
    for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
      send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
      for (auto &r : send_regions[d]) {
        x_low = 2 * r.x_low;
        y_low = r.y_low;
        z_low = r.z_low;
        x_high = 2 * r.x_high;
        y_high = r.y_high;
        z_high = r.z_high;
    
        if(x_low <= surface_x && surface_x < x_high && y_low <= surface_y && surface_y < y_high && z_low <= surface_z && surface_z < z_high){
          if(d == comm::DIM_X){
            auto it_from = exchange_surface_x[sect].find(surface_id);
            if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(surface_id);
          }else if(d == comm::DIM_Y) {
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

void ABVIModel::clear_exchange_ghost(){
  for (auto& vec : exchange_ghost) {
    vec.clear();
  }
}

void ABVIModel::clear_exchange_surface(const unsigned int next_sect){
  exchange_surface_x[next_sect].clear();
  exchange_surface_y[next_sect].clear();
  exchange_surface_z[next_sect].clear();
}

// std::array<std::vector<Lattice>, 7> ABVIModel::get_exghost(){
//   return exchange_ghost;
// }

// const comm::ColoredDomain* ABVIModel::get_pdomain(){
//   return p_domain;
// }

unsigned long ABVIModel::defectSize() {
  return 1; // todo add implementation
}