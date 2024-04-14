//
// Created by genshen on 2019/10/10.
//

#ifndef MISA_KMC_MODEL_ADAPTER_H
#define MISA_KMC_MODEL_ADAPTER_H

#include <comm/domain/colored_domain.h>
#include <type_define.h>
#include <utils/random/random.h>
#include <vector>
#include <array>
#include "lattice/lattice.h"
#include <unordered_set>
#include <comm/preset/sector_forwarding_region.h>

/**
 * \brief this is kmc model adapter
 * \tparam E the event type
 */
template <typename E> class ModelAdapter {
public:
  typedef comm::Region<comm::_type_lattice_size> lat_region;

  const comm::ColoredDomain *p_domain = nullptr;

  std::array<std::vector<ChangeLattice>, 7> exchange_ghost;

  // 8 个扇区的 surface 的按需发送区
  std::array<std::unordered_set<_type_lattice_id>, 8> exchange_surface_x;

  std::array<std::unordered_set<_type_lattice_id>, 8> exchange_surface_y;

  std::array<std::unordered_set<_type_lattice_id>, 8> exchange_surface_z;

  ModelAdapter();

  /**
   * \brief calculate rates in a region in the given region
   * \return the sum of all rates.
   */
  virtual _type_rate calcRates(const lat_region region) = 0;

  virtual _type_rate calcRatesGPU(const lat_region region, bool *sector_first, const lat_region region_next, int sect, int is_last_step) = 0;

  virtual void transFirstsectLocalToGpu(const lat_region region) = 0;

  virtual void clear_exchange_ghost() = 0;

  virtual void clear_exchange_surface(const unsigned int next_sect) = 0;

  // virtual _type_rate calcRatesGPU(const lat_region region, bool *sector_first) = 0;
  /**
   * \brief select an event from rates list in the given region
   * \param excepted_rate excepted rate
   * \param sum_rates the total rates
   */
  virtual E select(const lat_region region, const _type_rate excepted_rate, const _type_rate sum_rates) = 0;

  virtual unsigned long defectSize() = 0;

  /**
   * \brief perform the kmv event.
   */
  virtual void perform(const E e, const lat_region region, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) = 0;

  virtual void selectAndPerformOnGPU(const _type_rate rate, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) = 0;

  /**
   * \brief select an event and perform the event.
   * this is the wrapper function for \fn select() and \fn perform()
   * \param sum_rates the total rates
   */
  void selectAndPerform(const lat_region region, const _type_rate sum_rates, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id);

  void selectAndPerformGPU(int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id);

  /**
   * \brief reindex defect list in a specific region.
   * \param region the region
   */
  virtual void reindex(const lat_region region) = 0;

  /**
   * \brief generate a float random number between [0,1)
   * \return
   */
  inline double rand() {
    return r::random(); // todo random number
  }
};

template <typename E> ModelAdapter<E>::ModelAdapter() {}

template <class E> void ModelAdapter<E>::selectAndPerform(const lat_region region, const _type_rate sum_rates, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) {
  // [0, sum_rates)
  E e = select(region, rand() * sum_rates, sum_rates);
  // E e = select(rand() * sum_rates);
  // E e = select(region, rand() * 1.0, sum_rates);
  perform(e, region, rank, step, sect, sector_id, next_sector_id);
}

template <class E> void ModelAdapter<E>::selectAndPerformGPU(int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) {
  selectAndPerformOnGPU(rand(), rank, step, sect, sector_id, next_sector_id);
}

#endif // MISA_KMC_MODEL_ADAPTER_H
