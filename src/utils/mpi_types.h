//
// Created by genshen on 2019/10/26.
//

#ifndef MISA_KMC_MPI_TYPES_H
#define MISA_KMC_MPI_TYPES_H

#include <mpi.h>

namespace mpi_types {
  extern MPI_Datatype _mpi_type_lattice_data;
  extern MPI_Datatype _mpi_type_lattice_change;

  // called in initialization
  void setInterMPIType();

  // called when program finished
  void unsetInterMPIType();

  void setMPI_DataTypeLattice(MPI_Datatype *mpi_type_lat, MPI_Datatype *mpi_type_lat_change);
} // namespace mpi_types

#endif // MISA_KMC_MPI_TYPES_H
