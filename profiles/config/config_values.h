//
// Created by genshen on 2019/9/18.
//

#ifndef MISA_KMC_CONFIG_VALUES_H
#define MISA_KMC_CONFIG_VALUES_H

#include "lattice/lattice_types.h"
#include <comm/types_define.h>
#include <utils/bundle.h>
#include <vector>

namespace conf {
  enum CreateOption { None, Random, Pipe, Restart };
  struct Create {
    CreateOption create_option = CreateOption::None;
    ;
    // random
    int64_t va_count = 0;
    std::vector<LatticeTypes::lat_type> types;
    std::vector<int64_t> types_ratio;
    // pipe
    std::string pipe_input_box = "";
    // restart
    std::string restart_file = "";
  };

  struct Output {
    int64_t dump_interval = 0;
    std::string dump_file_path = "";
    int64_t logs_interval = 0;
    std::string logs_file = "";
    bool logs_to_file = false;
  };

  struct RandomSeeds {
    int64_t create_types;
    int64_t create_vacancy;
    int64_t event_selection;
    int64_t time_inc;
  };

  struct ConfigValues {
    friend std::ostream &operator<<(std::ostream &os, const ConfigValues &cv);

  public:
    Create create;
    Output output;
    RandomSeeds seeds;
    // box
    int64_t box_size[comm::DIMENSION_SIZE];

    double lattice_const = 0.0;

    double cutoff_radius = 0.0;

    // simulation
    double temperature = 0.0;
    double physics_time = 0.0;
    int64_t steps_limit = 0;
    double attempt_freq = 0.0;
    // isgenr in config file.
    bool is_def_gen = false;
    // dpasm1 in config file
    double dpa_ps = 0.0;

    void packData(kiwi::Bundle &bundle);

    void unpackData(kiwi::Bundle &bundle);
  };
} // namespace conf

#endif // MISA_KMC_CONFIG_VALUES_H
