//
// Created by runchu on 2019-09-24.
//

#include "simulation.h"
#include "abvi/kmc.h"
#include "pack/ghost_init_packer.h"
#include "utils/mpi_types.h"
#include "utils/simulation_domain.h"
#include <utils/mpi_utils.h>
#include <iostream>
// 进行一些领域分解和初始化操作，似乎与MPI通信有关。这段代码的目的是创建一个领域（Domain），并在其中进行一些设置。以下是对这段代码的简要解释：
// MPI_Comm new_comm;：声明一个MPI通信句柄，用于领域创建后的通信。
// comm::mpi_process pro = comm::mpi_process{...};：创建一个mpi_process结构，用于存储MPI进程的信息，包括本地进程的排名（rank）、所有排名的列表和通信句柄。
// int64_t phase_space_int64[comm::DIMENSION_SIZE];：声明一个数组phase_space_int64，用于存储phase_space数组的int64_t类型的副本。
// for (int i = 0; i < comm::DIMENSION_SIZE; i++) { ... }：将phase_space数组的每个元素转换为int64_t类型，并存储在phase_space_int64数组中。
// _p_domain = comm::ColoredDomain::Builder()...：使用comm::ColoredDomain的构建器模式创建一个领域（Domain）对象。
// 在这个过程中，设置了一些参数，包括通信句柄、相位空间、截断半径、晶格常数等。
// kiwi::mpiUtils::onGlobalCommChanged(new_comm);：更新全局MPI通信，以便将新的通信句柄应用于全局通信。这可能是在领域创建后更新通信配置的一部分。
// SimulationDomain::setSimDomain(kiwi::mpiUtils::global_process);：将模拟领域的信息设置为全局MPI进程的信息，这似乎是为了在后续模拟中使用这些信息。
// 这段代码的主要目的是初始化一个领域，并设置与领域相关的参数和通信句柄。如果你有关于特定部分的更具体的问题或需要更详细的解释，请提出相关问题。

void simulation::createDomain(const int64_t phase_space[comm::DIMENSION_SIZE], const double lattice_const,
                              const double cutoff_radius) {

  //进行区域分解
  //    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initializing GlobalDomain
  //    decomposition.\n");

  MPI_Comm new_comm;
  // 这里的local_process.own_rank，local_process.all_ranks， local_process.comm 在之前初始化过
  // kiwi::mpiUtils::local_process.comm = MPI_COMM_WORLD
  // 将 kiwi::mpiUtils类中的静态变量 local_process 复制给 comm::mpi_process 的 pro 这里的 local_process = global_process
  comm::mpi_process pro = comm::mpi_process{
      kiwi::mpiUtils::local_process.own_rank,
      kiwi::mpiUtils::local_process.all_ranks,
      kiwi::mpiUtils::local_process.comm,
  };
  // In domain creation, we set ghost size as (cut_lattice +1) to avoid
  // neighbor index overflowing. 在域创建中，我们将重影大小设置为（cut_lattice+1）以避免邻居索引溢出。
  int64_t phase_space_int64[comm::DIMENSION_SIZE];
  // 将 unsigned long 转成 int64 来存
  for (int i = 0; i < comm::DIMENSION_SIZE; i++) {
    phase_space_int64[i] = (int64_t)phase_space[i]; //就是box_size
  }

  // simulation 类中的成员变量 comm::ColoredDomain *_p_domain = nullptr;
  _p_domain = comm::ColoredDomain::Builder()
                  // _mpi_pro = pro; _p_comm = new_comm; 将_p_domain的Builder成员变量的两个数据初始化
                  .setComm(pro, &new_comm)
                  // _phase_space = phase_space_int64; 将_p_domain的Builder成员变量的_phase_space数据初始化
                  .setPhaseSpace(phase_space_int64)
                  // _cutoff_radius_factor = cutoff_radius / lattice_const 
                  // _ghost_size = static_cast<int>(ceil(_cutoff_radius_factor));
                  // 将_p_domain的Builder成员变量的 _cutoff_radius_factor 和 _ghost_size 数据初始化
                  // 即 _cutoff_radius_factor = 5.6 / 2.85532 = 1.96125
                  // 即 _ghost_size = 2
                  .setCutoffRadius(cutoff_radius / lattice_const)
                  // _lattice_const = lattice_const 将_p_domain的Builder成员变量的 _lattice_const 数据初始化
                  .setLatticeConst(lattice_const)
                  //            .setGhostSize(static_cast<int>(ceil(cutoff_radius/lattice_const)))
                  // 将_p_domain的父类Domain的成员变量 _grid_size 存储 三维拓扑的进程空间的三个维度的进程数量
                  // 将_p_domain的父类Domain的成员变量 _p_comm = 新创建的Cartesian拓扑通信器的指针。
                  // 将_p_domain的父类Domain的成员变量 _grid_coord = 当前进程在新创建的Cartesian拓扑通信器中的坐标(3维)
                  // 将_p_domain的父类Domain的成员变量 _rank_id_neighbours[3][2] = 每个维度上相邻进程的rank (在Cartesian拓扑通信器下) 
                  // _rank_id_neighbours[d][0] 表示低方向的邻近 在 x 轴为左邻近 d = 0,1,2 对应 x,y,z _rank_id_neighbours[d][1] 表示高方向的邻近 在 x 轴为右邻近
                  // 将_p_domain的父类Domain的成员变量 _meas_global_length[d] = _phase_space[d] * _lattice_const 即 box_size * 晶格常数存储在_meas_global_length中
                  // 将_p_domain的父类Domain的成员变量 _meas_global_region.low[d] = 0; // lower bounding is set to 0 by default.
                  // 将_p_domain的父类Domain的成员变量 _meas_global_region.high[d] = _meas_global_length[d];
                  // 将_p_domain的父类Domain的成员变量 _sub_box_lattice_size[d] = 当前进程在各个维度分到的 sub_box 的大小 
                  // 若 box 在 x 维度为 11, _grid_size 在 x 维度为 3 , 则 _grid_coord[0] = 0,1 的进程分到 4 _grid_coord[0] = 2 的进程分到 3
                  // 将_p_domain的父类Domain的成员变量 _lattice_size_ghost[d] = _ghost_size 即 _lattice_size_ghost[0,1,2] = 2
                  // 将_p_domain的父类Domain的成员变量 _ghost_extended_lattice_size[d] = _sub_box_lattice_size[d] + 2 * _lattice_size_ghost[d];
                  // 即 _ghost_extended_lattice_size[d] = _sub_box_lattice_size[d] + 4
                  // 将_p_domain的父类Domain的成员变量 _sub_box_lattice_region.low[d] = 存储的是当前进程拥有的子box的在整个box中的各个维度的起始索引。
                  // 将_p_domain的父类Domain的成员变量 _sub_box_lattice_region.high[d] = _sub_box_lattice_region.low[d] + _sub_box_lattice_size[d]; 存储的是当前进程拥有的子box的在整个box中的各个维度的终止索引。
                  // 将_p_domain的父类Domain的成员变量 _ghost_ext_lattice_region.low[d] = 当前进程拥有的子box的在整个box中各个轴的起始索引 - 2; 即扩大 2, 因为要包含ghost区域
                  // 将_p_domain的父类Domain的成员变量 _ghost_ext_lattice_region.high[d] = 当前进程拥有的子box的在整个box中各个轴的终止索引 + 2; 即扩大 2, 因为要包含ghost区域
                  // 将_p_domain的父类Domain的成员变量 _local_ghost_ext_lattice_region.low[d] = 0; 局部包含 ghost 区的sub_box区域
                  // 将_p_domain的父类Domain的成员变量 _local_ghost_ext_lattice_region.high[d] = _ghost_extended_lattice_size[d]; 局部包含 ghost 区的sub_box区域大小
                  // _local_ghost_ext_lattice_region.low[d] - _local_ghost_ext_lattice_region.high[d] 是当前进程的局部索引，加上ghost区后所负责区域维度的大小
                  // 将_p_domain的父类Domain的成员变量 _local_sub_box_lattice_region.low[d] = _lattice_size_ghost[d]; 局部不含 ghost 区的sub_box区域, 从 2 开始
                  // 将_p_domain的父类Domain的成员变量 _local_sub_box_lattice_region.high[d] = _lattice_size_ghost[d] + _sub_box_lattice_size[d];
                  // _local_sub_box_lattice_region.low[d] -- _local_sub_box_lattice_region.high[d] 是当前进程的局部索引，不加上ghost区后所负责区域维度的大小
                  // 将_p_domain的父类Domain的成员变量 _meas_sub_box_region.low[d] = 当前进程拥有的子box的在整个box中各个维度的起始索引 * 晶格常数
                  // 将_p_domain的父类Domain的成员变量 _meas_sub_box_region.high[d] = 当前进程拥有的子box的在整个box中的各个维度的终止索引 * 晶格常数
                  // 将_p_domain的父类Domain的成员变量 _meas_ghost_length[d] = 2 * 晶格常数 这个 2 就是 ghost 区各个维度的 size
                  // 将_p_domain的父类Domain的成员变量 _meas_ghost_ext_region.low[d] = _meas_sub_box_region.low[d] - 2 * 晶格常数
                  // 将_p_domain的父类Domain的成员变量 _meas_ghost_ext_region.high[d] = _meas_sub_box_region.high[d] + 2 * 晶格常数

                  // sector_lattice_size 一共有8个元素，每个元素都是一个包含三个int的数组 即 sector_lattice_size[8][3]
                  // 将_p_domain的成员变量 sector_lattice_size[0,2,4,6][0] = 当前进程在 x 轴维度分到的 sub_box 的大小 / 2;
                  // 将_p_domain的成员变量 sector_lattice_size[1,3,5,7][0] = 当前进程在 x 轴维度分到的 sub_box 的大小 / 2 - 当前进程在 x 轴维度分到的 sub_box 的大小 / 2;
                  // 将_p_domain的成员变量 sector_lattice_size[0,1,4,5][1] = 当前进程在 y 轴维度分到的 sub_box 的大小 / 2;
                  // 将_p_domain的成员变量 sector_lattice_size[2,3,6,7][1] = 当前进程在 y 轴维度分到的 sub_box 的大小 / 2 - 当前进程在 y 轴维度分到的 sub_box 的大小 / 2;
                  // 将_p_domain的成员变量 sector_lattice_size[0,1,2,3][2] = 当前进程在 z 轴维度分到的 sub_box 的大小 / 2;
                  // 将_p_domain的成员变量 sector_lattice_size[4,5,6,7][2] = 当前进程在 z 轴维度分到的 sub_box 的大小 / 2 - 当前进程在 z 轴维度分到的 sub_box 的大小 / 2;
                  // local_sector_region[8] 是一个一维数组，其中元素为 Region<int>
                  // local_sector_region[0] = Region( x_low = 0, y_low = 0, z_low = 0, x_high = sector_lattice_size[0][0], y_high = sector_lattice_size[0][1], z_high = sector_lattice_size[0][2])
                  // local_sector_region[1] = Region( x_low = 当前进程在 x 轴维度分到的 sub_box 的大小 / 2, y_low = 0, z_low = 0, x_high = 当前进程在 x 维度分到的大小, y_high = sector_lattice_size[0][1], z_high = sector_lattice_size[0][2])
                  // local_sector_region[2] = Region( x_low = 0, y_low = 当前进程在 y 轴维度分到的 sub_box 的大小 / 2, z_low = 0, x_high = sector_lattice_size[0][0], y_high = 当前进程在y维度分到的大小, z_high = sector_lattice_size[0][2])
                  // 其他local_sector_region同理
                  // 然后再将 local_sector_region[0 - 7].x_low -- z_high = local_sector_region[0 - 7].x_low -- z_high + 2; 
                  // local_ghost_ext_sector_region[8] 是一个一维数组，其中元素为 Region<int>
                  // local_ghost_ext_sector_region[0 - 7].x_low -- z_low = local_sector_region[0 - 7].x_low -- z_low - 2; 即为初始的 x_low
                  // local_ghost_ext_sector_region[0 - 7].x_high -- z_high = local_sector_region[0 - 7].x_high -- z_high + 2; 即为初始的 x_high + 4
                  // int local_split_coord[3]
                  // local_split_coord[0] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][0] + 2; 即 (当前进程在 x 轴维度分到的 sub_box 的大小 / 2) + 2;
                  // local_split_coord[1] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][1] + 2; 即 (当前进程在 y 轴维度分到的 sub_box 的大小 / 2) + 2;
                  // local_split_coord[2] = sector_lattice_size[X_LOW | Y_LOW | Z_LOW][2] + 2; 即 (当前进程在 z 轴维度分到的 sub_box 的大小 / 2) + 2;
                  .build(); 
  // 将 kiwi::mpiUtils 中的 global_process 设置为 Cartesian 拓扑通信器，并且也更新 global_process 的 all_ranks 和 own_rank 为当前进程在 Cartesian 拓扑通信器中的 rank 信息
  kiwi::mpiUtils::onGlobalCommChanged(new_comm); // set new domain.
  // 将 kiwi::mpiUtils 中的 global_process 复制给 SimulationDomain 结构体中的静态变量 kiwi_sim_pro 和 comm_sim_pro
  SimulationDomain::setSimDomain(kiwi::mpiUtils::global_process);

  //    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initialization done.\n");
}

void simulation::createLattice() {
  BoxBuilder builder{};
  // todo type conversion. 各种转换。
  // 每个进程设置自己的box_size为子box的晶格尺寸的维度
  // 将 builder 的 box_x, box_y, box_z 成员变量 = 当前进程在各个维度分到的 sub_box 的大小 
  // std::cout << " box_x is : " << _p_domain->sub_box_lattice_size[0] << " box_y is : " << _p_domain->sub_box_lattice_size[1] << " box_z is : " << _p_domain->sub_box_lattice_size[2] << std::endl;
  builder.setBoxSize(static_cast<_type_box_size>(_p_domain->sub_box_lattice_size[0]),
                     static_cast<_type_box_size>(_p_domain->sub_box_lattice_size[1]),
                     static_cast<_type_box_size>(_p_domain->sub_box_lattice_size[2]));
  // 设置该进程中的重影大小的x、y和z维度。 目前lattice_size_ghost 都为 2
  // 将 builder 的 ghost_x, ghost_y, ghost 成员变量 = 2 
  builder.setGhostSize(static_cast<_type_box_size>(_p_domain->lattice_size_ghost[0]),
                       static_cast<_type_box_size>(_p_domain->lattice_size_ghost[1]),
                       static_cast<_type_box_size>(_p_domain->lattice_size_ghost[2]));
  // 为生成器设置全局模拟框的晶格大小， 记录整个大的box_size (80,80,80)
  // 将 builder 的 global_box_x, global_box_y, global_box_z 成员变量 = 整个 box 在三个维度的大小
  builder.setGlobalLatSize(_p_domain->phase_space[0], _p_domain->phase_space[1], _p_domain->phase_space[2]);
  // _sub_box_lattice_region.low[0] = 当前进程的x轴坐标 * (box_size[x] / 笛卡尔域中的x轴进程数) + min(当前进程的x轴坐标, box_size[x] % 笛卡尔域中的x轴进程数)
  // 为生成器设置全局模拟框的起始晶格坐标
  // 将 builder 的 global_base_x, global_base_y, global_base_z 成员变量 = 当前进程拥有的子box的在整个box中的各个维度的起始索引
  builder.setGlobalBaseLat(_p_domain->sub_box_lattice_region.x_low, _p_domain->sub_box_lattice_region.y_low,
                           _p_domain->sub_box_lattice_region.z_low);
  // create empty lattice list(lattice types are not specified) and defects 创建空晶格列表（未指定晶格类型）和缺陷
  // list.
  // 就是新 new 了一个 box, 然后以上面初始化过的 box_x, box_y, box_z, ghost_x, ghost_y, ghost, global_box_x, global_box_y, global_box_z, global_base_x, global_base_y, global_base_z 作为参数初始化一个 LatListMeta 对象
  // 然后用 LatListMeta 对象去初始化 box 的 lattice_list, va_list, itl_list 。 lattice_list = new PeriodLatticeList(lattice_meta);
  // LatListMeta 中有成员变量 size_x, size_y, size_z; box_x, box_y, box_z; g_box_x, g_box_y, g_box_z; g_base_x, g_base_y, g_base_z; ghost_x, ghost_y, ghost_z; _max_id; local_base_id
  // size_x = (2 * 当前进程在x维度分到的 sub_box 的大小  + 2 * 2 * ghost_x), size_y = (当前进程在y维度分到的 sub_box 的大小  + 2 * ghost_y), size_z = (当前进程在z维度分到的 sub_box 的大小  + 2 * ghost_z) 
  // 这里的 size 实际上就是当前进程在各个维度分到的 sub_box 的大小加上两边的 ghost区
  // box_x = (2 * 当前进程在x维度分到的 sub_box 的大小 ), box_y = (当前进程在y维度分到的 sub_box 的大小 ), box_z = (当前进程在z维度分到的 sub_box 的大小 )
  // g_box_x = (2 * box 在x维度的大小), g_box_y = (box 在y维度的大小), g_box_z = (box 在z维度的大小)
  // g_base_x = (2 * 当前进程拥有的子box的在整个box中的x维度的起始索引), g_base_y = (当前进程拥有的子box的在整个box中的y维度的起始索引), g_base_z = (当前进程拥有的子box的在整个box中的z维度的起始索引)
  // ghost_x = (2 * ghost_x), ghost_y = (ghost_y), ghost_z = (ghost_z)
  // _max_id = size_z * size_y * size_x - 1
  box = builder.build(); // todo delete pointer
}

void simulation::prepareForStart() {
  // initialize ghost area for all sectors. 初始化所有子域的重影区域。
  GhostInitPacker packer{_p_domain, box->lattice_list};
  // 在计算之前进程将自己的一部分数据发给邻居进程，这些数据作为邻居进程的ghost区。
  // 其他就是在更新各个进程的ghost区，使得ghost区内的数据与邻居进程真正的数据保持一致 
  // comm::neiSendReceive(&packer, SimulationDomain::comm_sim_pro, mpi_types::_mpi_type_lattice_data,
  //                      _p_domain->rank_id_neighbours);
  const int DoubleSideForwardingTag2 = 0x100;

  const MPI_Datatype data_type = packer.getMPI_DataType();

  int num_send[comm::DIMENSION_SIZE][2];
  int num_receive[comm::DIMENSION_SIZE][2];
  Lattice *send_buff[2];
  Lattice *receive_buff[2];

  MPI_Status status;
  MPI_Status send_statuses[comm::DIMENSION_SIZE][2];
  MPI_Status recv_statuses[comm::DIMENSION_SIZE][2];
  MPI_Request send_requests[comm::DIMENSION_SIZE][2];
  MPI_Request recv_requests[comm::DIMENSION_SIZE][2];

  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    int dimension = d;
    // 循环两次
    for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
      num_send[dimension][direction] = packer.sendLength2(dimension, direction);
      assert(num_send[dimension][direction] >= 0); // todo remove
      // send_buff[0]是一个指针，指向了一个新分配的lattice数组，该数组的大小为num_send[0][0]
      send_buff[direction] = new Lattice[num_send[dimension][direction]];
      packer.onSend2(send_buff[direction], num_send[dimension][direction], dimension, direction);
    }
    for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
      int &numsend = num_send[dimension][direction];
      int numrecv = 0;

      MPI_Isend(send_buff[direction], numsend, data_type,
                packer.p_domain->rank_id_neighbours[dimension][direction], DoubleSideForwardingTag2,
                SimulationDomain::comm_sim_pro.comm, &send_requests[dimension][direction]);
      // test the status of neighbor process.
      MPI_Probe(packer.p_domain->rank_id_neighbours[dimension][(direction + 1) % 2],
                DoubleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &status);
      // test the data length to be received.
      MPI_Get_count(&status, data_type, &numrecv);
      // initialize receive buffer via receiving size.
      // the receiving length is get bt MPI_Probe from its neighbour process.
      assert(numrecv >= 0);  // todo remove
      receive_buff[direction] = new Lattice[numrecv];
      num_receive[dimension][direction] = numrecv;
      MPI_Irecv(receive_buff[direction], numrecv, data_type,
                packer.p_domain->rank_id_neighbours[dimension][(direction + 1) % 2], DoubleSideForwardingTag2,
                SimulationDomain::comm_sim_pro.comm, &recv_requests[dimension][direction]);
    }
    for (int direction = comm::DIR_LOWER; direction <= comm::DIR_HIGHER; direction++) {
      MPI_Wait(&send_requests[dimension][direction], &send_statuses[dimension][direction]);
      MPI_Wait(&recv_requests[dimension][direction], &recv_statuses[dimension][direction]);
      packer.onReceive2(receive_buff[direction], num_receive[dimension][direction], dimension, direction);
      // release buffer
      delete[] send_buff[direction];
      delete[] receive_buff[direction];
    }
  }
  // all finished
  packer.onFinish();
}

void simulation::get_neighbour_local_sub_box() {
  const int DoubleSideForwardingTag2 = 0x1100;
  MPI_Status send_statuses[comm::DIMENSION_SIZE];
  MPI_Status recv_statuses[comm::DIMENSION_SIZE];
  MPI_Request send_requests[comm::DIMENSION_SIZE];
  MPI_Request recv_requests[comm::DIMENSION_SIZE];
  _type_lattice_size send_buff[comm::DIMENSION_SIZE];
  // int receive_buff[3];
  for (int dimension = 0; dimension < comm::DIMENSION_SIZE; dimension++) {
    send_buff[dimension] = _p_domain->sub_box_lattice_size[dimension];
    // 循环两次
    MPI_Isend(send_buff + dimension, 1, MPI_UNSIGNED_LONG_LONG,
              _p_domain->rank_id_neighbours[dimension][1], DoubleSideForwardingTag2,
              SimulationDomain::comm_sim_pro.comm, &send_requests[dimension]);

    MPI_Irecv(_p_domain->neighbour_local_sub_box + dimension, 1, MPI_UNSIGNED_LONG_LONG,
              _p_domain->rank_id_neighbours[dimension][0], DoubleSideForwardingTag2,
              SimulationDomain::comm_sim_pro.comm, &recv_requests[dimension]);
          
  }
  for (int dimension = 0; dimension < comm::DIMENSION_SIZE; dimension++) {
    MPI_Wait(&send_requests[dimension], &send_statuses[dimension]);
    MPI_Wait(&recv_requests[dimension], &recv_statuses[dimension]);
  }
}