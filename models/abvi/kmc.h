//
// Created by genshen on 2018-12-12.
//

#ifndef MISA_KMC_KMC_H
#define MISA_KMC_KMC_H

#include <vector>
#include <array>
#include "lattice/lattice.h"
#include "box.h"
#include "event.h"
#include "plugin/event_listener.h"
#include "type_define.h"
#include <models/model_adapter.h>
#include <comm/domain/colored_domain.h>
#include <unordered_map>
#include "abvi/defect/vac_hash.h"
#include <comm/preset/sector_forwarding_region.h>
#include "utils/simulation_domain.h"
#include <comm/comm.hpp>
#include <cmath>

// #include "../../src/algorithms/sl/sublattice.h"

// #include "../gpu/gpu_simulate.h"
/*!
 * \brief the model routine of KMC simulation, including rate calculation, event selecting
 * and execution implementation.
 */

class ABVIModel : public ModelAdapter<event::SelectedEvent> {

public:
  /**
   * \brief initialize kmc with simulation box.
   */
  explicit ABVIModel(Box *box, double v, double T);

  /**
   * \brief calculate the transition rates of each defect
   *  in a region at current process
   *
   * The calculating methods depends on the lattice type(single atom, vacancy and dumbbell).
   * Different types of lattice have different methods or formulas to calculate the rate.
   * see the implementation for more details.
   *
   * \param region a region, this function will calculate the transition rates in this region
   * \note the x lattice size is not doubled in \param region parameter
   * After this step, the rate of every transition direction of each lattice will be set.
   * \return return the sum of rates of all KMC events in this region,
   *  including dumbbell transition and vacancy transition and defect generation.
   */
  _type_rate calcRates(const lat_region region) override; // todo


  _type_rate calcRatesGPU(const lat_region region, bool *sector_first, const lat_region region_next, int sect, int is_last_step) override;
  // _type_rate calcRatesGPU(const lat_region region, bool *sector_first) override; // todo
  // _type_rate calcRates_GPU(const lat_region region);

  void transFirstsectLocalToGpu(const lat_region region) override;

  /*
   * \brief select an event randomly from rates list in a given region.  从给定区域的机率列表中随机选择一个事件。
   *
   * \param excepted_rate which equals to total rate* random number between 0-1. 等于总机率*0-1之间的随机数。
   * \param total_rates the sum rates 机率之和
   * \note the x lattice size is not doubled in \param region parameter \param region参数中的x晶格大小没有加倍
   * \return the selected event.  返回所选事件。
   */
  event::SelectedEvent select(const lat_region region, const _type_rate excepted_rate,
                              const _type_rate sum_rates) override;

  /**
   * \brief perform the selected KMC event.
   *
   */
  void perform(const event::SelectedEvent event, const lat_region region, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) override;

  void selectAndPerformOnGPU(const _type_rate rate, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) override;

  void reindex(const lat_region region) override;

  /**
   * \brief set kmc event listener.
   * \param p_listener pointer to the event listener.
   */
  void setEventListener(EventListener *p_listener);

  void setColoredDomain(comm::ColoredDomain *_p_domain);

  void addExchange_ghost(Lattice lat_to, const unsigned int sector_id);

  void addExchange_surface(_type_lattice_id surface_id);

  void clear_exchange_ghost() override;

  void clear_exchange_surface(const unsigned int  next_sect) override;

  unsigned long defectSize() override;

  int add_test_commu(const unsigned int sector_id);

  int add_test_surface_commu(const unsigned int next_sector_id);

protected:
  double time = 0;

public:
  Box *box = nullptr; // todo init box pointer
  /**
   * \brief attempt frequency.
   */
  const double v;
  /**
   * \brief temperature
   */
  const double T; 

  /**
   * \brief pointer to event listener.
   * event callback function will be called when executing a kmc event.
   */
  EventListener *p_event_listener = nullptr;

  /**
   * \brief it returns the rate of defect generation.
   * \return
   */
  _type_rate defectGenRate();

  // const comm::ColoredDomain *p_domain = nullptr;


};

#endif // MISA_KMC_KMC_H
