//
// Created by genshen on 2019/9/17.
//

#include "pkmc.h"
#include <sys/time.h>
#include <iostream>

int main(int argc, char *argv[]) {
  struct timeval time_start, time_end;
  double timeuse;
  gettimeofday(&time_start,NULL);
  (new PKMC())->run(argc, argv);
  gettimeofday(&time_end,NULL);
  timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e6;
  // //if(SimulationDomain::comm_sim_pro.own_rank % 1000 == 0 ){
  std::cout << " misa-akmc total time is : "<< timeuse << " s " << std::endl;
  //}
  return 0;
}
