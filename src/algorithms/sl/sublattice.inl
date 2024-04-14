//
// Created by genshen on 2019/10/12.
//

#ifndef MISA_KMC_SUB_LATTICE_INC
#define MISA_KMC_SUB_LATTICE_INC

#include "comm_dirs.h"
#include "utils/random/random.h"
#include "utils/simulation_domain.h"
#include "pack/sim_sync_packer.h"
#include <cmath>
#include <comm/comm.hpp>
#include <comm/preset/sector_forwarding_region.h>
#include <logs/logs.h>
#include <sys/time.h>
#include <iostream>
#include <mpi.h>
#include "../gpu/gpu_simulate.h"

//#include <sys/sysinfo.h>

template <class PKg, class PKs, class Ins, typename E>
void SubLattice::startTimeLoop(Ins pk_inst, ModelAdapter<E> *p_model, EventHooks *p_event_hooks) {
  const int64_t time_steps = ceil(time_limit / T); // T 是写死的3e-7，time_limit是config.yaml文件中的physics_time
  double time_total_start, time_total_end;
  double time_calculate_start, time_calculate_end;
  double time_commu_start, time_commu_end;
  double max_time_commu, time_calculate_total, time_commu_total;
  time_commu_total = 0;
  time_calculate_total = 0;
  double vac_time_total = 0.0;
  // struct sysinfo s_info;
  //int error0;
  // gpu_prepare(const long size_x, const long size_y, const long size_z, const double v, const double T)
  bool sector_first;
  bool sector_first2;
  int is_last_step = 0;
  //int ii = 0;
  int itera1 = 0;
  int itera2 = 0;
  int max_itera1 = 0;
  int max_itera2 = 0;
  // p_event_hooks->onStepFinished(0);
  if(SimulationDomain::comm_sim_pro.own_rank == 0) {
    kiwi::logs::v(" ", " rank id is : {} x dimension divided is : {} y dimension divided is : {} z dimension divided is : {}.\n", 
    SimulationDomain::comm_sim_pro.own_rank, 
    p_model->p_domain->sub_box_lattice_size[0],
    p_model->p_domain->sub_box_lattice_size[1],
    p_model->p_domain->sub_box_lattice_size[2]);
  }

  //time_total_start = MPI_Wtime();
  //p_model->transFirstsectLocalToGpu(p_domain->local_sector_region[(*sec_meta.sector_itl).id]);
  //time_total_end = MPI_Wtime() - time_total_start;
  //if(SimulationDomain::comm_sim_pro.own_rank == 0) time_calculate_total += time_total_end;
  // 扇区的执行顺序是 0, 7, 2, 5, 3, 4, 1, 6
  MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
  time_total_start = MPI_Wtime();
  for (int64_t step = 0; step < time_steps; step++) { // time steps loop
    // transFirstsectToGpu(p_model, (*sec_meta.sector_itl).id);
    if(step == time_steps - 1) is_last_step = 1;
    for (int sect = 0; sect < SECTORS_NUM; sect++) {      // sector loop SECTORS_NUM = 8 依据同步子域算法，在每个进程的区域内再进行子域的划分(一个8个子域)
    //step_threshold_time = (step + 1) * 3.0e-7 - 演化时间 evolution_time 开始的时候是0
      //test_gpu();
      const double step_threshold_time = static_cast<double>(step + 1) * T - sec_meta.sector_itl->evolution_time; 
      // if(SimulationDomain::comm_sim_pro.own_rank == 0 && step > 52 && step < 56){
      //   kiwi::logs::v(" ", "rank id is : {} step is : {} sect is : {} step_threshold_time is : {} sec_meta.sector_itl->evolution_time is : {} T is : {}.\n", SimulationDomain::comm_sim_pro.own_rank, step, sect, step_threshold_time, sec_meta.sector_itl->evolution_time, T);
      // }
      double sector_time = 0.0;
      // p_model->reindex(p_domain->local_sector_region[(*sec_meta.sector_itl).id]); // reindex0 defects lists
      p_model->clear_exchange_ghost();
      p_model->clear_exchange_surface((*sec_meta.sector_itl).id);
      sector_first = true;
      sector_first2 = true;
      //ii = 0;
      //const type_sector next_sector = sec_meta.sector_itl.next();
      //if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " step is : {} sect is : {} starting !!!.\n", step, sect);
      // MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
      while (sector_time < step_threshold_time) { // note: step_threshold_time may be less then 0.0
        // 各个方向的迁移机率，并将其累加
        time_calculate_start = MPI_Wtime();
        const double total_rates = calcRatesWrapper(p_model, (*sec_meta.sector_itl).id, &sector_first, sect, is_last_step);
        time_calculate_end = MPI_Wtime() - time_calculate_start;
        if(SimulationDomain::comm_sim_pro.own_rank == 0) vac_time_total += time_calculate_end;
        // if(SimulationDomain::comm_sim_pro.own_rank == 0 && sector_first2) {
        //   kiwi::logs::v(" ", " step is : {} sect is : {} total_rates is : {}.\n", step, sect, total_rates);
        //   sector_first2 = false;
        // }
        // if(total_rates != 0 && sector_first2){
        //   kiwi::logs::v(" ", "rank is : {} step is : {} sect is : {} step_threshold_time is : {} sector_time is : {} total_rates is : {} i is : {} .\n", SimulationDomain::comm_sim_pro.own_rank, step, sect, step_threshold_time, sector_time, total_rates, i);
        //   // std::cout << " rank is : " << SimulationDomain::comm_sim_pro.own_rank << " step is : " << step << " sect is : " << sect << " step_threshold_time is : " << step_threshold_time << " sector_time is :" << sector_time << " total_rates is : " << total_rates << std::endl;
        //   sector_first2 = false;
        // }
        //ii++;
        // if(total_rates != 0 && sector_first2){
        //   std::cout << " rank is : " << SimulationDomain::comm_sim_pro.own_rank << " step is : " << step << " sect is : " << sect << " step_threshold_time is : " << step_threshold_time << " sector_time is :" << sector_time << " total_rates is : " << total_rates << std::endl;
        // }
        // std::numeric_limits<_type_rate>::epsilon() 是用于获取类型 _type_rate 的最小可表示正数的函数。
        // 在C++中，std::numeric_limits 是一个模板类，它提供了与特定数据类型相关的各种信息，包括数据类型的最小值、最大值和精度等。
        // epsilon() 函数用于获取该类型能够表示的最小正数值，即最接近零的正数。
        // 所以，std::abs(total_rates) < std::numeric_limits<_type_rate>::epsilon() 表达式的含义是：检查 total_rates 的绝对值是否小于类型 _type_rate 的最小正数值。
        // 如果这个条件成立，那么 total_rates 被认为足够接近零，即它在数值上被视为零或接近零。
        // 这种比较通常用于处理浮点数的数值精度问题。如果 total_rates 的绝对值小于 _type_rate 类型的最小正数值，
        // 那么可以认为 total_rates 在精度上可以被视为零，这在一些数值计算和比较中是有用的。
        if (total_rates == 0.0 || std::abs(total_rates) < std::numeric_limits<_type_rate>::epsilon()) {
          // If there is no defect, use synchronous parallel kMC
          // algorithm. Because there is no kMC event, just increase
          // time. 如果没有缺陷，则使用同步并行kMC算法。因为没有kMC事件，只是增加时间。这样便直接跳出了循环
          sector_time = step_threshold_time;
        } else {
          // 有KMC事件
          selectPerformWrapper(p_model, total_rates, (*sec_meta.sector_itl).id, SimulationDomain::comm_sim_pro.own_rank, step, sect);
          // time_inc_dis can produces random numbers in a range [0, 1)
          // time_inc_rng 是config.yaml中给出的随机数种子
          const double rand = time_inc_dis(time_inc_rng);
          const double delta_t = -std::log(rand) / total_rates;
          sector_time += delta_t;
        }
        // todo: time comparing, nearest principle based on predicting
        // next delta t.  时间比较，基于预测下一个Δt的最接近原理。
      }
      // kiwi::logs::v(" ", "rank is : {} step is : {} sect is : {} i is : {} .\n", SimulationDomain::comm_sim_pro.own_rank, step, sect, i);
      // gettimeofday(&time_end,NULL);
      // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e6;
      // std::cout << " rank is : " << SimulationDomain::comm_sim_pro.own_rank << " calculate time is : " << timeuse <<  " s " << std::endl;
      // MPI_Reduce(&timeuse, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, SimulationDomain::comm_sim_pro.comm);
      if(!sector_first) sector_final();
      // 不断的发生事件，当sector_time >= step_threshold_time后结束。
      // kmc simulation on this sector is finished
      // we store evolution time of current sector when the sector loop
      // finishes.  该扇区的kmc模拟已完成，当扇区循环结束时，我们存储当前扇区的演变时间。
      sec_meta.sector_itl->evolution_time += sector_time;
      // communicate ghost area of current process to sync simulation
      // regions of neighbor process. 
      // 通信当前进程的重影区域以同步相邻进程的模拟区域。
      // gettimeofday(&time_end,NULL);
      // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0e3 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e3;
      // std::cout << "rank is : " << SimulationDomain::comm_sim_pro.own_rank << "step is :  " << step << "sect is : " << sect << "time is : " << timeuse << std::endl;
      // gettimeofday(&time_start,NULL);
      // int ghost_commu_size = p_model->add_test_commu((*sec_meta.sector_itl).id);
      // int surface_commu_size = p_model->add_test_surface_commu(next_sector.id);
      //MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
      //calculateNextRegionAdv(sect);
      time_commu_start = MPI_Wtime();
      syncSimRegions<PKs>(pk_inst, p_model->exchange_ghost, p_model->exchange_surface_x, p_model->exchange_surface_y, p_model->exchange_surface_z);
      syncNextSectorGhostRegions<PKg>(pk_inst, p_model->exchange_surface_x, p_model->exchange_surface_y, p_model->exchange_surface_z); // communicate ghost area data of next sector in
                                                // current process. 通信当前进程中下一个扇区的重影区域数据。

      //MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
      time_commu_end = MPI_Wtime() - time_commu_start;
      MPI_Reduce(&time_commu_end, &max_time_commu, 1, MPI_DOUBLE, MPI_MAX, 0, SimulationDomain::comm_sim_pro.comm);
      if(SimulationDomain::comm_sim_pro.own_rank == 0) time_commu_total += max_time_commu;
      // gettimeofday(&time_end,NULL);
      // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0e3 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e3;
      // std::cout << "rank is : " << SimulationDomain::comm_sim_pro.own_rank << "step is :  " << step << "sect is : " << sect << "com_time is : " << timeuse << std::endl;
      ++sec_meta.sector_itl;                    // update sector id. 更新扇区id。
      nextSector();                             // some post operations after moved to next sector. 在移动到下一个子域后的一些后期操作。目前是空的
      // gettimeofday(&time_end2,NULL);
      // timeuse2 = (time_end2.tv_sec - time_start2.tv_sec) * 1.0 + (double)(time_end2.tv_usec - time_start2.tv_usec) / 1.0e6;
      // std::cout << " rank is : " << SimulationDomain::comm_sim_pro.own_rank << " communicate time is : " << timeuse2 <<  " s " << std::endl;
      // error0 = sysinfo(&s_info);
      // std::cout << "rank is : " << SimulationDomain::comm_sim_pro.own_rank << "error0 :  " << error0 << "totalram : " << s_info.totalram / (1024.0 * 1024 * 1024) << "freeram : " << s_info.freeram / (1024.0 * 1024 * 1024) << std::endl;

    }
    p_event_hooks->onStepFinished(step);
  }
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) time_calculate_total += time_total_end;
  p_event_hooks->onAllDone(); //空的函数
  gpu_final();
  MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
  time_total_end = MPI_Wtime() - time_total_start;
  //print_kernel_time();
  // MPI_Reduce(&itera1, &max_itera1, 1, MPI_INT, MPI_MAX, 0, SimulationDomain::comm_sim_pro.comm);
  // MPI_Reduce(&itera2, &max_itera2, 1, MPI_INT, MPI_MAX, 0, SimulationDomain::comm_sim_pro.comm);
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", "  max_itera1  1 is : {} . \n", max_itera1);
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", "  max_itera2  2 is : {} . \n", max_itera2);
  if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " vac time is : {} total time is : {} s. communicate time is : {} s. \n", vac_time_total, time_total_end, time_commu_total);
}

template <typename E> double SubLattice::calcRatesWrapper(ModelAdapter<E> *p_model, const type_sector_id sector_id, bool *sector_first, int sect, int is_last_step) {
   const type_sector next_sector = sec_meta.sector_itl.next();
   return p_model->calcRatesGPU(p_domain->local_sector_region[sector_id], sector_first, p_domain->local_sector_region[next_sector.id], sect, is_last_step);
   //return p_model->calcRates(p_domain->local_sector_region[sector_id]);
}

template <typename E>
void SubLattice::selectPerformWrapper(ModelAdapter<E> *p_model, const _type_rate total_rates,
                                      const type_sector_id sector_id, int rank, _type_lattice_count step, int sect) {
  // 进程当前的子域id 和 各个方向的迁移机率之和作为参数
  const type_sector next_sector = sec_meta.sector_itl.next();
  p_model->selectAndPerformGPU(rank, step, sect, sector_id, next_sector.id);
  //p_model->selectAndPerform(p_domain->local_sector_region[sector_id], total_rates, rank, step, sect, sector_id, next_sector.id);
}

// PKs：打包器（Packer）的类型，用于打包数据以进行通信。
// Ins：通信实例，用于管理通信的配置和状态。
// 获取当前模拟区域的部分信息：
// 获取当前模拟区域的扇区（cur_sector）。
// 获取通信维度的顺序，这里通信顺序是Z、Y、X。
// 计算发送和接收区域：
// 为每个通信维度（X、Y、Z）计算发送区域（send_regions）和接收区域（recv_regions）。
// 使用comm::fwCommSectorRecvRegion和comm::fwCommSectorSendRegion函数计算，根据当前扇区、维度、幽灵尺寸等参数计算。
// 将计算的区域信息存储在send_regions和recv_regions数组中。
// 计算通信方向和排名：
// 计算发送方向和接收方向的数组，分别为send_dirs和recv_dirs。
// 计算发送和接收方向的排名，分别为ranks_send和ranks_recv。
// 创建打包器对象：
// 创建用于打包数据的PKs类型的打包器对象packer。
// 通信顺序是从Z到Y再到X。
// 执行单向前向通信：
// 调用comm::singleSideForwardComm函数进行通信，实际的通信操作在这里完成。
// 函数参数包括打包器对象packer、通信进程、数据类型、发送区域、接收区域、发送排名和接收排名等信息。
// 释放资源：
// 在每个通信维度上，发送和接收数据后，释放发送缓冲区和接收缓冲区的内存。
// 完成通信：
// 调用打包器对象的onFinish方法，完成通信操作。
// 总之，syncSimRegions函数实现了在不同通信维度上的单向前向通信，通过计算发送和接收区域，将数据从一个区域传输到另一个区域，并使用打包器管理数据的打包和解包。
// 这种通信操作有助于在并行计算中同步模拟区域的数据。



// class PackerInstance {
// public:
//   explicit PackerInstance(LatticesList *lattice_list) : lattice_list(lattice_list) {}

//   SimSyncPacker newSimCommPacker() { return SimSyncPacker{lattice_list}; }

//   GhostSyncPacker newGhostCommPacker() { return GhostSyncPacker{lattice_list}; }

// private:
//   // reference or pointers to create packer
//   LatticesList *lattice_list;
// };




template <class PKs, class Ins> void SubLattice::syncSimRegions(Ins &pk_inst, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, 
                                                                std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                                                std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y, 
                                                                std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z) {
  // Note that, in this communication, the sending regions are equals to 注意，在该通信中，发送区域等于在每个维度同步重影时使用的接收区域。
  // the receiving regions used in synchronizing ghost at each dimension.
  // Besides, the communication order of dimensions is Z,Y,X; 此外，维度的传播顺序为Z、Y、X；
  // but in synchronizing ghost communication, the order is X,Y,Z. 但是在同步重影通信中，顺序是X、Y、Z。
  // todo the regions can be static. 到目前为止，这些区域可以是静态的。

  const type_sector cur_sector = (*sec_meta.sector_itl);
  const type_sector next_sector = sec_meta.sector_itl.next();

  // const int dims[3] = {0, 1, 2};
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z}; 

  // typedef std::array<std::vector<comm::Region<int>>, 3> type_comm_lat_regions;
  // 声明了一个数组，数组大小为3，数组中的各个元素为vector，vector中存储的是comm::Region<int>类型的数据。

  type_comm_lat_regions send_regions; // send regions in each dimension
  type_comm_lat_regions recv_regions; // receive regions in each dimension.

  // for (int d = 0; d < 3; d++) {
  // 根据sector_id确定哪些部分的数据需要发送(确定范围)
  // 根据sector_id确定哪些部分需要接收数据(确定范围)
  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    // cur_sector.id是当前扇区的id
    send_regions[d] = comm::fwCommSectorRecvRegion(cur_sector.id, dims[d], p_domain->lattice_size_ghost,
                                                   p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
    recv_regions[d] = comm::fwCommSectorSendRegion(cur_sector.id, dims[d], p_domain->lattice_size_ghost,
                                                   p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
  }

  // pack data 
  // call ssfdCommRecvDirs for send dirs and call ssfdCommSendDirs for receive
  // 调用ssfdCommRecvDirs获取发送指令，调用ssfdCommeSendDirs获取接收指令
  // dirs.

  // 实际上就是把0~7的数转化为二级制的形式，如sector_id为5，则return一个{1，0，1}
  const std::array<unsigned int, comm::DIMENSION_SIZE> send_dirs = ssfdCommRecvDirs(cur_sector.id);
  // 实际上就是把7-(0~7)的数转化为二级制的形式，如sector_id为5，则return一个{0，1，0}
  const std::array<unsigned int, comm::DIMENSION_SIZE> recv_dirs = ssfdCommSendDirs(cur_sector.id);

  // 根据sector_id判断要发送到哪些邻居进程，将3个邻居进程的rank_id存于ranks_send数组中
  const std::array<unsigned int, comm::DIMENSION_SIZE> ranks_send = {
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_X][send_dirs[0]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Y][send_dirs[1]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Z][send_dirs[2]]),
  };
  // 根据sector_id判断要接收哪些邻居进程，将3个邻居进程的rank_id存于ranks_recv数组中
  const std::array<unsigned int, comm::DIMENSION_SIZE> ranks_recv = {
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_X][recv_dirs[0]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Y][recv_dirs[1]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Z][recv_dirs[2]]),
  };

  // 初始化一个packer,貌似是用于打包数据进程进程通信的？
  // int m = 0;
  // for (auto& vec : exchange_ghost) {
  //   m += vec.size();
  // }
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " ghost commu size is : {} .\n", m);
  // for (size_t i = 0; i < exchange_ghost.size(); ++i) {
  //   // 对每个元素，遍历std::vector<ChangeLattice>中的每个ChangeLattice对象
  //   for (size_t j = 0; j < exchange_ghost[i].size(); ++j) {
  //     // 访问当前元素
  //     ChangeLattice& currentLattice = exchange_ghost[i][j];
  //     kiwi::logs::v(" ", " x is : {} y is : {} z is : {}.\n", currentLattice.x, currentLattice.y, currentLattice.z);
  //   }
  // }
  // kiwi::logs::v("", "1111111111111111111111111111111111111111111111111111 .\n");
  // for(int i = 0; i < 8; i++){
  //   for (const auto& lattice : exchange_ghost[i]) {
  //     kiwi::logs::v(" ", " type is : {} x is : {} y is : {} z is : {}.\n", lattice.type._type, lattice.x, lattice.y, lattice.z);
  //     // buffer[len++] = lattice;
  //   }
  // }
  const int SingleSideForwardingTag2 = 0x101;
  const int SendCountTag2 = 0x110;

  SimSyncPacker packer = pk_inst.newSimCommPacker();
  // note here, communication order is from Z to Y, and X. 注意，这里的通信顺序是从Z到Y，以及X。
  // pack_date_type = Lattice      pack_region_type = int
  // static inline MPI_Datatype getMPI_DataType() { return mpi_types::_mpi_type_lattice_data; }
  // comm::singleSideForwardComm<typename PKs::pack_date_type, typename PKs::pack_region_type, true>(
  //     &packer, SimulationDomain::comm_sim_pro, packer.getMPI_DataType(), send_regions, recv_regions, ranks_send,
  //     ranks_recv);
  const MPI_Datatype data_type = packer.getMPI_DataTypeChange();
  int num_send[comm::DIMENSION_SIZE];
  int num_receive[comm::DIMENSION_SIZE];
  ChangeLattice *send_buff;
  ChangeLattice *receive_buff;

  MPI_Status status;

  MPI_Status send_statuses[comm::DIMENSION_SIZE];
  MPI_Status recv_statuses[comm::DIMENSION_SIZE];
  MPI_Request send_requests[comm::DIMENSION_SIZE];
  MPI_Request recv_requests[comm::DIMENSION_SIZE];

  ChangeLattice send_count;
  for (int d = 2; d > -1; d--) {
    int dim_id = d;
    // 里的通信顺序是从Z到Y，以及X。
    // if (F) {
    //     dim_id = DIMENSION_SIZE - 1 - d;
    // }
    // dim_id = comm::DIMENSION_SIZE - 1 - d;
    // if()
    // if(ranks_send[dim_id] != SimulationDomain::comm_sim_pro.own_rank){

    // }
    num_send[dim_id] = packer.sendLength2(dim_id, exchange_ghost, exchange_surface_x[next_sector.id], &send_count);
    assert(num_send[dim_id] >= 0); // todo remove
    //if(num_send[dim_id] > 0){
    send_buff = new ChangeLattice[num_send[dim_id]];
    // 将这次要发送的数据复制到send_buff
    packer.onSend2(send_buff, exchange_ghost, dim_id, exchange_surface_x[next_sector.id], &send_count);

    // send and received data.
    int &numsend = num_send[dim_id];

    // 将send_buff中的数据单向发送到send_ranks[dim_id]进程中
    MPI_Isend(send_buff, numsend, data_type, ranks_send[dim_id], SingleSideForwardingTag2,
              SimulationDomain::comm_sim_pro.comm, &send_requests[dim_id]);
    // MPI_Barrier(SimulationDomain::comm_sim_pro.comm);
    // kiwi::logs::v(" ", "1111111111111111111111111111 .\n" );
    int numrecv = 0;
    // int flag;
    // test the status of neighbor process. 测试邻居进程的状态
    //MPI_Iprobe(ranks_recv[dim_id], SingleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &flag, &status);
    MPI_Probe(ranks_recv[dim_id], SingleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &status);
    //kiwi::logs::v(" ", "333333333333333333333333333333333 .\n" );
    //assert(status.MPI_SOURCE == ranks_recv[dim_id]);
    //assert(status.MPI_TAG == SingleSideForwardingTag2);
    // test the data length to be received. 获取要接收的数据的量的大小
    //if(flag){
    MPI_Get_count(&status, data_type, &numrecv);
    // initialize receive buffer via receiving size.
    // the receiving length is get bt MPI_Probe from its neighbour process.
    assert(numrecv >= 0);  // todo remove
    //kiwi::logs::v(" ", "2222222222222222222222222222222 .\n");
    //if(numrecv > 0){
    // 从recv_ranks[dim_id]进程接收数据到receive_buff中
    receive_buff = new ChangeLattice[numrecv];
    num_receive[dim_id] = numrecv;
    MPI_Irecv(receive_buff, numrecv, data_type, ranks_recv[dim_id],
              SingleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &recv_requests[dim_id]);

    MPI_Wait(&send_requests[dim_id], &send_statuses[dim_id]);
    MPI_Wait(&recv_requests[dim_id], &recv_statuses[dim_id]);
    //将接收到的数据存到_lattices中，也就是更新全局的lattices信息
    packer.onReceive2(receive_buff, recv_regions[dim_id], num_receive[dim_id], dim_id, exchange_ghost, 
                      p_domain->sub_box_lattice_size, p_domain->neighbour_local_sub_box, cur_sector.id, 
                      exchange_surface_x, exchange_surface_y, exchange_surface_z, p_domain);
    // release buffer
    //if(num_send[dim_id] > 0) delete[] send_buff;
    delete[] send_buff;
    delete[] receive_buff;
    //if(numrecv > 0) delete[] receive_buff;
  
  }
    // all finished
    packer.onFinish();
}

// 这是一个用于同步下一个扇区（sector）的幽灵（ghost）区域数据的C++函数。
// 函数的目的是在模拟中处理幽灵区域的数据通信，以确保相邻扇区之间的数据保持同步。以下是函数的主要参数和步骤：
// PKg packer：用于打包和解包数据的打包器对象。
// Ins &pk_inst：打包器对象的实例，用于创建通信打包器。
// type_sector next_sector：下一个扇区的信息，用于确定通信目标。
// const int dims[comm::DIMENSION_SIZE]：包含通信维度的数组，例如X、Y和Z维度。
// type_comm_lat_regions send_regions：包含发送区域信息的数组，用于不同维度。
// type_comm_lat_regions recv_regions：包含接收区域信息的数组，用于不同维度。
// std::array<unsigned int, comm::DIMENSION_SIZE> send_dirs：包含发送方向的数组，用于不同维度。
// std::array<unsigned int, comm::DIMENSION_SIZE> recv_dirs：包含接收方向的数组，用于不同维度。
// std::array<unsigned int, comm::DIMENSION_SIZE> ranks_send：包含发送排名的数组，用于指定数据发送到哪个进程。
// std::array<unsigned int, comm::DIMENSION_SIZE> ranks_recv：包含接收排名的数组，用于指定从哪个进程接收数据。
// 函数的主要步骤如下：
// 根据下一个扇区的信息，计算在每个通信维度上的发送和接收区域信息。这些区域信息存储在send_regions和recv_regions数组中。
// 根据通信维度和扇区的信息，计算发送和接收方向的数组send_dirs和recv_dirs。
// 根据发送和接收方向，确定发送和接收的进程排名，并将它们存储在ranks_send和ranks_recv数组中。
// 创建一个通信打包器对象packer，用于打包和解包数据。这个打包器对象是通过pk_inst.newGhostCommPacker()创建的。
// 调用comm::singleSideForwardComm函数，执行单向前向通信，将数据从发送区域传输到接收区域。在这个过程中，使用打包器对象来打包和解包数据。
// 总之，这个函数用于在不同通信维度上同步下一个扇区的幽灵区域数据，确保模拟中的相邻扇区之间的数据保持同步。这对于并行模拟中的数据一致性非常重要。

template <class PKg, class Ins> void SubLattice::syncNextSectorGhostRegions(Ins &pk_inst, 
                                                                            std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                                                            std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y, 
                                                                            std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z) {
  const type_sector cur_sector = (*sec_meta.sector_itl);
  const type_sector next_sector = sec_meta.sector_itl.next();
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  // todo the regions can be static.
  type_comm_lat_regions send_regions; // send regions in each dimension
  type_comm_lat_regions recv_regions; // receive regions in each dimension.
  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    send_regions[d] = comm::fwCommSectorSendRegion(next_sector.id, dims[d], p_domain->lattice_size_ghost,
                                                   p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
    recv_regions[d] = comm::fwCommSectorRecvRegion(next_sector.id, dims[d], p_domain->lattice_size_ghost,
                                                   p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
  }

  // pack data
  
  // 实际上就是把7-(0~7)的数转化为二级制的形式，如sector_id为5，则return一个{0，1，0}
  const std::array<unsigned int, comm::DIMENSION_SIZE> send_dirs = ssfdCommSendDirs(next_sector.id);
  // 实际上就是把0~7的数转化为二级制的形式，如sector_id为5，则return一个{1，0，1}
  const std::array<unsigned int, comm::DIMENSION_SIZE> recv_dirs = ssfdCommRecvDirs(next_sector.id);
  const std::array<unsigned int, comm::DIMENSION_SIZE> ranks_send = {
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_X][send_dirs[0]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Y][send_dirs[1]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Z][send_dirs[2]]),
  };
  const std::array<unsigned int, comm::DIMENSION_SIZE> ranks_recv = {
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_X][recv_dirs[0]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Y][recv_dirs[1]]),
      static_cast<unsigned int>(p_domain->rank_id_neighbours[comm::DIM_Z][recv_dirs[2]]),
  };


  //if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " surface commu size is : {} .\n", exchange_surface_x[next_sector.id].size() + exchange_surface_y[next_sector.id].size() + exchange_surface_z[next_sector.id].size());

  // do data pack
  PKg packer = pk_inst.newGhostCommPacker();
  // comm::singleSideForwardComm(&packer, SimulationDomain::comm_sim_pro, packer.getMPI_DataType(), send_regions,
  //                             recv_regions, ranks_send, ranks_recv);

  const int SingleSideForwardingTag2 = 0x101;
  // const MPI_Datatype data_type = packer.getMPI_DataTypeChange();
  const MPI_Datatype data_type = packer.getMPI_DataType();
  int num_send[comm::DIMENSION_SIZE];
  int num_receive[comm::DIMENSION_SIZE];
  Lattice *send_buff;
  Lattice *receive_buff;

  MPI_Status status;
  MPI_Status send_statuses[comm::DIMENSION_SIZE];
  MPI_Status recv_statuses[comm::DIMENSION_SIZE];
  MPI_Request send_requests[comm::DIMENSION_SIZE];
  MPI_Request recv_requests[comm::DIMENSION_SIZE];

  for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
    int dim_id = d;
    // prepare data
    // note: the direction in sendLength, onSend and onReceive are not used (ignored).
    //num_send[dim_id] = packer.sendLength2(dim_id, exchange_surface_x[next_sector.id], exchange_surface_y[next_sector.id], exchange_surface_z[next_sector.id]);
    num_send[dim_id] = packer.sendLength(send_regions[dim_id], dim_id, comm::DIR_LOWER);
    assert(num_send[dim_id] >= 0); // todo remove
    send_buff = new Lattice[num_send[dim_id]];
    //packer.onSend2(send_buff, num_send[dim_id], dim_id, exchange_surface_x[next_sector.id], exchange_surface_y[next_sector.id], exchange_surface_z[next_sector.id]);
    packer.onSend(send_buff, send_regions[dim_id], num_send[dim_id], dim_id, comm::DIR_LOWER);
    // send and received data.
    int &numsend = num_send[dim_id];
    int numrecv = 0;
    MPI_Isend(send_buff, numsend, data_type, ranks_send[dim_id], SingleSideForwardingTag2,
              SimulationDomain::comm_sim_pro.comm, &send_requests[dim_id]);
    // test the status of neighbor process.
    MPI_Probe(ranks_recv[dim_id], SingleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &status);
    // test the data length to be received.
    MPI_Get_count(&status, data_type, &numrecv);
    // initialize receive buffer via receiving size.
    // the receiving length is get bt MPI_Probe from its neighbour process.
    assert(numrecv >= 0);  // todo remove
    receive_buff = new Lattice[numrecv];
    num_receive[dim_id] = numrecv;
    MPI_Irecv(receive_buff, numrecv, data_type, ranks_recv[dim_id],
              SingleSideForwardingTag2, SimulationDomain::comm_sim_pro.comm, &recv_requests[dim_id]);

    // data received.
    MPI_Wait(&send_requests[dim_id], &send_statuses[dim_id]);
    MPI_Wait(&recv_requests[dim_id], &recv_statuses[dim_id]);
    // packer.onReceive2(receive_buff, num_receive[dim_id], dim_id, exchange_surface_x, exchange_surface_y, exchange_surface_z,
    //                   p_domain->sub_box_lattice_size, p_domain->neighbour_local_sub_box, cur_sector.id, p_domain);
    packer.onReceive(receive_buff, recv_regions[dim_id], num_receive[dim_id], dim_id, comm::DIR_LOWER);
    // release buffer
    delete[] send_buff;
    delete[] receive_buff;
  }

  // all finished
  packer.onFinish();
}

#endif // MISA_KMC_SUB_LATTICE_INC
