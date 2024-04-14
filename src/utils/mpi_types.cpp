//
// Created by genshen on 2019/10/26.
//

#include "mpi_types.h"
#include "lattice/lattice.h"
#include <cassert>

MPI_Datatype mpi_types::_mpi_type_lattice_data;
MPI_Datatype mpi_types::_mpi_type_lattice_change;

void mpi_types::setInterMPIType() { 
  setMPI_DataTypeLattice(&_mpi_type_lattice_data, &_mpi_type_lattice_change); 
}

void mpi_types::unsetInterMPIType() { 
  MPI_Type_free(&_mpi_type_lattice_data); 
  MPI_Type_free(&_mpi_type_lattice_change); 
}

// const int entries = 2;：定义了数据结构中的成员数量，这里是2，因为有两个成员。
// MPI_Datatype array_of_types[] = {MPI_INT, MPI_UNSIGNED_LONG};：定义了一个数组，其中包含了数据结构中每个成员的 MPI 数据类型。一个成员是 MPI_INT 类型，另一个是 MPI_UNSIGNED_LONG 类型。
// int array_of_block_len[] = {1, 1};：定义了一个数组，其中包含了数据结构中每个成员的块长度。在这里，每个成员都有长度1，表示每个成员都占据一个连续的内存块。
// Lattice lat_dummy;：创建一个 Lattice 类型的临时对象 lat_dummy，以便获取对象中各个成员的内存地址。
// 使用条件编译和 MPI 版本检查：
// MPI_Get_address 函数用于获取对象 lat_dummy 及其成员 type 和 id 的内存地址。
// 使用 assert 断言检查获取的地址，以确保它们在内存中的顺序和关系正确。
// 在满足 MPI 版本要求的情况下，使用 MPI_Type_create_struct 创建一个自定义的 MPI 数据类型，该数据类型描述了 Lattice 数据结构的布局。
// 使用 MPI_Type_commit 来提交并完成 MPI 数据类型的创建。
// 总之，这个函数的主要目的是在 MPI 中创建一个自定义的数据类型，该数据类型用于描述 Lattice 数据结构的布局，以便在 MPI 进程之间传递 Lattice 类型的数据。
// 这种做法可以确保数据在不同进程之间正确地传递和解释。此外，函数中的条件编译和 MPI 版本检查确保了代码在满足特定 MPI 版本要求的情况下才会执行。

void mpi_types::setMPI_DataTypeLattice(MPI_Datatype *mpi_type_lat, MPI_Datatype *mpi_type_lat_change) {
  const int entries = 2;
  const int total = 4;
  // Lattice::type is int type and Lattice::id is unsigned long type.
  MPI_Datatype array_of_types[] = {MPI_INT, MPI_UNSIGNED_LONG_LONG};
  // type, x, y, z
  MPI_Datatype array_of_types_ghost[] = {MPI_INT, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG};
  // 1 int type and 1 unsigned long type.
  int array_of_block_len[] = {1, 1};
  int array_of_block_len_ghost[] = {1, 1, 1, 1};

  Lattice lat_dummy;
  ChangeLattice changelattice;
// MPI_Aint 是 MPI（消息传递接口）标准中定义的数据类型，用于表示内存地址的整数类型。它是 MPI 中的一种跨平台的整数类型，通常用于描述内存地址和位移量。
// MPI_Get_address 获取内存中某个位置的地址
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
  MPI_Aint addr_base, addr_type, addr_id;
  MPI_Aint addr_base_change, addr_type_change, addr_x, addr_y, addr_z;
  MPI_Get_address(&lat_dummy, &addr_base);
  MPI_Get_address(&lat_dummy.type, &addr_type);
  MPI_Get_address(&lat_dummy.id, &addr_id);

  MPI_Get_address(&changelattice, &addr_base_change);
  MPI_Get_address(&changelattice.type, &addr_type_change);
  MPI_Get_address(&changelattice.x, &addr_x);
  MPI_Get_address(&changelattice.y, &addr_y);
  MPI_Get_address(&changelattice.z, &addr_z);

  assert(addr_type < addr_id);
  assert(addr_base <= addr_type);

  assert(addr_base_change <= addr_type_change);
  assert(addr_type_change < addr_x);
  assert(addr_x < addr_y);
  assert(addr_y < addr_z);

  MPI_Aint array_of_dis[] = {addr_type - addr_base, addr_id - addr_base};

  MPI_Aint array_of_dis_change[] = {addr_type_change - addr_base_change, addr_x - addr_base_change, addr_y - addr_base_change, addr_z - addr_base_change};

  MPI_Type_create_struct(entries, array_of_block_len, array_of_dis, array_of_types, mpi_type_lat);
  MPI_Type_commit(mpi_type_lat);
  MPI_Type_create_struct(total, array_of_block_len_ghost, array_of_dis_change, array_of_types_ghost, mpi_type_lat_change);
  MPI_Type_commit(mpi_type_lat_change);
#else
  static_assert(false, "MPI version not matching (MPI 2 or above is needed)");
#endif
}
