//
// Created by genshen on 2019/9/18.
//

#include "config_parsing.h"
#include "lattice_types_string.h"
#include <comm/types_define.h>
#include <fstream>
#include <iostream>
#include <random>
#include <utils/random/random.h>

ConfigParsing *ConfigParsing::m_pInstance = nullptr;

ConfigParsing::ConfigParsing() : configValues() {}

ConfigParsing *ConfigParsing::newInstance(const std::string config_file) {
  if (m_pInstance == nullptr) {
    m_pInstance = new ConfigParsing(); // todo delete
    m_pInstance->parse(config_file);
  }
  return m_pInstance;
}

ConfigParsing *ConfigParsing::getInstance() {
  if (m_pInstance == nullptr) {
    m_pInstance = new ConfigParsing();
  }
  return m_pInstance; // make sure there is a configure instance.
}

bool ConfigParsing::configureCheck() { return false; }

void ConfigParsing::putConfigData(kiwi::Bundle &bundle) { configValues.packData(bundle); }

void ConfigParsing::getConfigData(kiwi::Bundle &bundle) { configValues.unpackData(bundle); }

// 打开指定路径的配置文件（config_file），并创建一个输入文件流（std::ifstream）来读取文件内容。
// 检查文件是否成功打开。如果无法打开文件，会设置错误消息，指示无法访问配置文件，并返回。
// 使用YAML库加载配置文件内容到YAML::Node对象中，该对象表示了配置文件的层次结构。
// 关闭文件流，不再需要访问文件。
// 解析配置文件的不同部分，包括：
// 解析名为"box"的部分，这似乎是关于盒子的配置信息。
// 解析名为"simulation"的部分，这似乎是关于模拟的配置信息。
// 解析名为"create"的部分，这似乎是关于创建的配置信息。
// 解析名为"seeds"的部分，这似乎是关于种子的配置信息。
// 解析名为"output"的部分，如果存在的话，这是有关输出的配置信息。
// 对于每个部分，首先检查是否存在该部分，如果不存在，则设置错误消息并返回。
// 对于每个部分，调用相应的解析函数（如parseBox，parseSim，parseCreate等）来解析该部分的配置信息。如果解析失败，则设置错误消息并返回。
// 最后，如果配置文件中存在"output"部分（这是可选的），则也尝试解析该部分的配置信息。如果解析失败，也会设置错误消息，但不会返回错误，因为"output"部分是可选的。
// 总之，这个函数的主要目的是解析指定的配置文件，将不同部分的配置信息加载到相应的数据结构中，以供程序后续使用。如果解析过程中出现问题，它会设置错误消息并返回，以便在调用代码中处理错误。
// 这是一个常见的配置文件解析函数，通常用于配置文件驱动的应用程序。

void ConfigParsing::parse(const std::string config_file) {
  std::ifstream ifs(config_file);
  if (!ifs.good()) { // todo fixme important, if file not exist, this branch
                     // would not enter.
    setError("can not access the configure file: " + config_file);
    return;
  }
  const YAML::Node config = YAML::Load(ifs);
  ifs.close();

  // parse box
  const YAML::Node box = config["box"];
  if (!box) {
    setError("box in config is not specified.");
    return;
  }
  if (!parseBox(box)) {
    return;
  }

  // parse simulation
  const YAML::Node sim = config["simulation"];
  if (!sim) {
    setError("simulation in config is not specified.");
    return;
  }
  if (!parseSim(sim)) {
    return;
  }

  // parse create.
  const YAML::Node create = config["create"];
  if (!sim) {
    setError("create in config is not specified.");
    return;
  }
  if (!parseCreate(create)) {
    return;
  }

  // parse random
  const YAML::Node seeds = config["seeds"];
  if (!sim) {
    setError("seeds in config is not specified.");
    return;
  }
  if (!parseSeeds(seeds)) {
    return;
  }

  // parse output
  const YAML::Node output = config["output"];
  if (output) { // it does not matter if output term in config file is not
                // specified.
    if (!parseOutput(output)) {
      return;
    }
  }
}

bool ConfigParsing::parseBox(const YAML::Node &yaml_box) {
  if (!yaml_box["size"].IsSequence() || yaml_box["size"].size() != comm::DIMENSION_SIZE) {
    setError("box size must be 3.");
    return false;
  } else {
    for (int i = 0; i < comm::DIMENSION_SIZE; i++) {
      configValues.box_size[i] = yaml_box["size"][i].as<int64_t>(0);
    }
  }

  if (!yaml_box["lattice_const"]) {
    setError("lattice const must be set in config.");
    return false;
  } else {
    configValues.lattice_const = yaml_box["lattice_const"].as<double>(0.0);
  }

  if (!yaml_box["cutoff_radius"]) {
    setError("cutoff radius must be set in config.");
    return false;
  } else {
    configValues.cutoff_radius = yaml_box["cutoff_radius"].as<double>(0.0);
  }
  return true;
}

bool ConfigParsing::parseSim(const YAML::Node &yaml_sim) {
  configValues.temperature = yaml_sim["temperature"].as<double>(0.0);
  configValues.physics_time = yaml_sim["physics_time"].as<double>(0.0);
  configValues.steps_limit = yaml_sim["steps_limit"].as<int64_t>(0);
  configValues.is_def_gen = yaml_sim["isgenr"].as<bool>(false);
  configValues.dpa_ps = yaml_sim["dpasm1"].as<double>(0.0);
  configValues.attempt_freq = yaml_sim["attempt_freq"].as<double>(0.0);
  return true;
}

bool ConfigParsing::parseCreate(const YAML::Node &yaml_create) {
  const YAML::Node yaml_random = yaml_create["random"];
  if (yaml_random) {
    const bool random_enabled = yaml_random["enable"].as<bool>(false);
    if (random_enabled) {
      configValues.create.create_option = conf::CreateOption::Random;
      configValues.create.va_count = yaml_random["va_count"].as<int64_t>(0);
      YAML::Node alloy = yaml_random["alloy"];
      if (alloy.IsMap()) {
        for (YAML::const_iterator it = alloy.begin(); it != alloy.end(); ++it) {
          configValues.create.types.emplace_back(lat::LatTypes(it->first.as<std::string>())); // .first是lattice.types.h文件中定义的 Fe = 0x0001；等作为map里面的key
          configValues.create.types_ratio.emplace_back(it->second.as<int64_t>(0)); // .second是config.yaml文件中定义的 Fe = 9598；等作为map里面的value,是合金内各个元素的比率
        }
      } else {
        setError("alloy must be a map in config.");
        return false;
      }
      return true;
    }
  }

  const YAML::Node yaml_pipe = yaml_create["pipe"];
  if (yaml_pipe) {
    const bool pipe_enabled = yaml_pipe["enable"].as<bool>(false);
    if (pipe_enabled) {
      configValues.create.create_option = conf::CreateOption::Pipe;
      configValues.create.pipe_input_box = yaml_pipe["input_box"].as<std::string>("");
      YAML::Node alloy = yaml_random["alloy"];
      if (alloy.IsMap()) {
        for (YAML::const_iterator it = alloy.begin(); it != alloy.end(); ++it) {
          configValues.create.types.emplace_back(lat::LatTypes(it->first.as<std::string>()));
          configValues.create.types_ratio.emplace_back(it->second.as<int64_t>(0));
        }
      } else {
        setError("alloy must be a map in config.");
        return false;
      }
      return true;
    }
  }

  const YAML::Node yaml_restart = yaml_create["restart"];
  if (yaml_restart) {
    const bool restart_enabled = yaml_restart["enable"].as<bool>(false);
    if (restart_enabled) {
      configValues.create.create_option = conf::CreateOption::Restart;
      configValues.create.restart_file = yaml_restart["file_path"].as<std::string>("");
      return true;
    }
  }

  if (configValues.create.create_option == conf::CreateOption::None) {
    setError("no create option in config.");
    return false;
  }
  return true;
}

bool ConfigParsing::parseSeeds(const YAML::Node &yaml_seeds) {
  configValues.seeds.create_types = yaml_seeds["create_types"].as<int64_t>(r::seed_auto);
  configValues.seeds.create_vacancy = yaml_seeds["create_vacancy"].as<int64_t>(r::seed_auto);
  configValues.seeds.event_selection = yaml_seeds["event_selection"].as<int64_t>(r::seed_auto);
  configValues.seeds.time_inc = yaml_seeds["time_inc"].as<int64_t>(r::seed_auto);
  return true;
}

bool ConfigParsing::parseOutput(const YAML::Node &output) {
  const YAML::Node yaml_dump = output["dump"];
  const YAML::Node yaml_logs = output["logs"];
  if (yaml_dump) {
    configValues.output.dump_interval = yaml_dump["interval"].as<int64_t>(0);
    configValues.output.dump_file_path = yaml_dump["file_path"].as<std::string>("");
  }
  if (yaml_logs) {
    configValues.output.logs_interval = yaml_logs["interval"].as<int64_t>(0);
    configValues.output.logs_file = yaml_logs["logs_file"].as<std::string>("");
  }
  // if it is empty, log to /dev/console
  configValues.output.logs_to_file = !(configValues.output.logs_file.empty());
  return true;
}
