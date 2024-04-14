//
// Created by runchu on 2019-09-24.
//

#ifndef MISA_KMC_SIMULATION_H
#define MISA_KMC_SIMULATION_H

#include <comm/domain/colored_domain.h>
#include <cstring>
#include <io/io_writer.h>
#include <mpi.h>
#include "utils/simulation_domain.h"
#include "abvi/box.h"
#include "algorithms/sl/sublattice.h"
#include "lattice/period_lattice_list.h"
#include "models/model_adapter.h"
#include "pack/ghost_sync_packer.h"
#include "pack/packer_instance.h"
#include "pack/sim_sync_packer.h"
#include <iostream>

class simulation {
public:
  simulation() = default;

  // virtual ~simulation();

  /**
   * Denote N as the count of all processors.
   * {@memberof domainDecomposition} will divide the simulation box into N
   * parts, we call each part as a sub-box. And each sub-box will bind to a
   * processor.
   * @param phase_space the lattice count in each dimension.
   * @param lattice_const lattice constance
   * @param cutoff_radius cutoff radius factor = cutoff/lattice_const
   */
  void createDomain(const int64_t phase_space[comm::DIMENSION_SIZE], const double lattice_const,
                    const double cutoff_radius);

  /**
   * \brief allocate memory for lattice list with ghost regions, and set ids
   * for each lattice site. \note the lattice types are not set here.
   */
  void createLattice();

  /**
   * Initialize ghost area for starting simulation.
   */
  void prepareForStart();

  void get_neighbour_local_sub_box();

  /**
   * \brief perform simulation using a specific model.
   * \tparam E type of event in kmc model.
   * \param p_model pointer of kmc model.
   * \param event_hooks pointer of event hooks or callbacks for handing each
   * algorithm event. \param seed_time_inc random number generation seed for
   * kmc time increasing. \param time_limit the max simulation time.
   */
  template <typename E>
  void simulate(ModelAdapter<E> *p_model, EventHooks *p_event_hooks, const double seed_time_inc,
                const double time_limit);

public:
  Box *box = nullptr;
  comm::ColoredDomain *_p_domain = nullptr;
};

template <typename E>
void simulation::simulate(ModelAdapter<E> *p_model, EventHooks *p_event_hooks, const double seed_time_inc,
                          const double time_limit) {
  SubLattice sl(_p_domain, seed_time_inc, time_limit,
                3.0e-7); // todo calculate T
  // if(SimulationDomain::comm_sim_pro.own_rank == 0){
  //   box->lattice_list->_lattices[2][2][4].type._type = LatticeTypes::V;
  //   box->va_list->mp.insert(std::make_pair(box->lattice_list->_lattices[2][2][4].id, Vacancy{}));
  // }
  PackerInstance pk_ins(box->lattice_list);
  // double time1, timeuse;
  // time1 = MPI_Wtime();
  sl.startTimeLoop<GhostSyncPacker, SimSyncPacker, PackerInstance>(pk_ins, p_model, p_event_hooks);
  // timeuse = MPI_Wtime() - time1;
  // std::cout << " simulate tiem is : " << timeuse << std::endl;
}

#endif // MISA_KMC_SIMULATION_H
