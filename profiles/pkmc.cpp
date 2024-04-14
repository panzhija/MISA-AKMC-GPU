//
// Created by genshen on 2019/9/18.
//

#include <abvi/kmc.h>
#include <args.hpp>
#include <iostream>
#include <lattice/lattices_list.h>
#include <logs/logs.h>
#include <utils/mpi_types.h>
#include <utils/mpi_utils.h>
#include <iostream>
#include "building_config.h"
#include "config/config_parsing.h"
#include "config/lattice_types_string.h"
#include "creation.h"
#include "device.h"
#include "m_event_hook.h"
#include "m_event_listener.h"
#include "pkmc.h"
#include "profile_config.h"
#include "seed_relocate.h"
#include "../gpu/gpu_simulate.h"

// 这个函数的主要功能是解析命令行参数，并根据解析结果执行不同的操作。让我解释一下函数的主要功能：
// 首先，它创建了一个名为parser的命令行参数解析器，使用了一个开源的C++库（args），用于帮助解析命令行参数。
// 然后，它定义了一些命令行参数的选项，包括：
// -h或--help：用于显示帮助菜单。
// -c或--conf：用于指定一个配置文件的路径。
// -v或--version：用于显示程序的版本号。
// 接下来，它使用parser.ParseCLI(argc, (const char *const *)argv)来解析传递给函数的命令行参数，其中argc是命令行参数数量，argv是命令行参数数组。
// 接着，它使用try-catch块来捕获可能发生的异常，包括：
// args::Help：如果用户请求帮助菜单，则打印帮助菜单并返回false，表示解析失败。
// args::ParseError：如果发生命令行参数解析错误，则打印错误消息和帮助菜单，并返回false，表示解析失败。
// args::ValidationError：如果发生命令行参数验证错误，则打印错误消息和帮助菜单，并返回false，表示解析失败。
// 如果解析成功，并且传递了-c或--conf参数，则将配置文件路径存储在configFilePath变量中，并返回true表示解析成功。
// 如果传递了-v或--version参数，则打印程序的版本信息和构建时间，并返回false表示解析成功。
// 最后，如果没有传递任何参数或解析失败，它会打印帮助菜单并返回false表示解析失败。
// 这个函数的主要目的是解析命令行参数，并根据参数执行不同的操作，例如显示帮助信息、显示版本信息或设置配置文件路径。如果解析成功，它返回true，否则返回false。

bool PKMC::beforeCreate(int argc, char **argv) {
  // parser arguments
  // see https://github.com/Taywee/args for using args.
  args::ArgumentParser parser("This is parallel misa-kmc program.", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<std::string> conf(parser, "conf", "The configure file", {'c', "conf"});
  args::Flag version(parser, "version", "show version number", {'v', "version"});
  try {
    parser.ParseCLI(argc, (const char *const *)argv);
  } catch (args::Help) {
    std::cout << parser;
    return false;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return false;
  } catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return false;
  }

  if (conf) {
    configFilePath = args::get(conf);
    return true;
  }

  if (version) {
    std::cout << "MISA KMC version " << KMC_VERSION_STRING << std::endl;
    std::cout << "Build time: " << __TIME__ << " " << __DATE__ << "." << std::endl;
    return false;
  }
  // if no args, print usage.
  std::cerr << parser;
  return false;
}

// 创建一个指向ConfigParsing类对象的指针p_config，该对象用于解析和存储配置信息。
// 如果当前进程是主进程（own_rank等于MASTER_PROCESSOR），则执行以下操作：
// 打印一条日志消息，表示MPI环境已经初始化。
// 创建一个ConfigParsing对象，通常是用来解析配置文件的。函数调用ConfigParsing::newInstance(configFilePath)，其中configFilePath是配置文件的路径。
// 检查配置对象是否有错误（hasError），如果有错误，则打印错误消息，然后调用this->abort(2)，可能是用于退出程序的函数。
// 如果当前进程不是主进程（不是主进程的工作节点），则简单地创建一个空的ConfigParsing对象。
// 使用配置对象（p_config）来同步配置数据到其他处理器，这可能涉及将配置信息从主进程传播到其他工作节点。
// 如果定义了编译时的宏KMC_DEBUG_MODE，则执行以下操作：
// 如果当前进程是主进程，将配置信息打印到标准输出（控制台）。
// 根据配置决定日志的设置：
// 如果程序输出到控制台并且控制台是真正的终端（没有I/O重定向），则启用彩色日志。
// 如果配置指定将日志输出到文件，则将日志输出到指定的文件。
// 总之，这个函数的主要任务是初始化配置对象，根据配置设置日志选项，以及在需要时将配置信息同步到其他MPI处理器。

void PKMC::onCreate() {
  // struct timeval time_start, time_end;
  // double timeuse;
  // gettimeofday(&time_start,NULL);
  ConfigParsing *p_config;
  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {  //own_rank 是当前进程在MPI_COMM_WORLD中的序号，all_ranks 是MPI_COMM_WORLD域中的所有进程数
    kiwi::logs::s("env", "mpi env is initialized.\n");
  }
  // initial config Obj, then read and resolve config file.
  if(kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) kiwi::logs::v(" ", " read config starting !!!.\n");
  p_config = ConfigParsing::newInstance(configFilePath); // todo config file from argv.
  if (p_config->hasError) {
    kiwi::logs::e("config", "{0}\n", p_config->errorMessage);
    this->abort(2);
  }
  //} else {
    // just initial a empty config Obj.
    //p_config = ConfigParsing::getInstance();
  //}
  // sync config data to other processors from master processor.
  if(kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) kiwi::logs::v(" ", " read config ending !!!.\n");
#ifdef KMC_DEBUG_MODE
  // print configure.
  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
    std::cout << p_config->configValues;
  }
#endif

  // prepare logs.
  if (istty() && !(p_config->configValues.output.logs_to_file)) {
    // set colorful log if we output to console and it is a real tty(no io
    // redirection).
    kiwi::logs::setColorFul(true);
  } else if (p_config->configValues.output.logs_to_file) {
    // 这里不执行
    kiwi::logs::setLogFile(p_config->configValues.output.logs_file);
  }
  // gettimeofday(&time_end,NULL);
  // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0e3 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e3;
  // std::cout << " rank : " << kiwi::mpiUtils::global_process.own_rank << " onCreate time : " << timeuse << " ms " << std::endl;
}

// mpi_types::setInterMPIType();：在 MPI 中创建一个自定义的数据类型，该数据类型用于描述 Lattice 数据结构的布局，以便在 MPI 进程之间传递 Lattice 类型的数据。
// sim = new simulation();：它创建了一个simulation类的实例，并将其分配给指针sim。这表明PKMC类包含了一个类型为simulation的成员变量。
// conf::ConfigValues config_v = ConfigParsing::getInstance()->configValues;：这一行从某个ConfigParsing类的单例实例中检索配置值，并将其存储在变量config_v中。很可能这些配置值在整个函数中都被使用。
// r::seedRelocate(...);：它似乎是重新定位用于随机数生成的种子值。具体的细节取决于r的实现和传递给它的参数。
// sim->createDomain(...);：这一行调用了sim对象上的createDomain方法，并传递了来自config_v的值作为参数。该方法很可能用于初始化simulation类中的一些与域相关的数据。
// sim->createLattice();：这个方法调用在simulation对象中初始化了一个晶格，这似乎是模拟中的一个重要部分。
// 基于config_v.create.create_option的switch语句：这个switch语句检查配置中的create_option的值，根据其值执行不同的操作：
// 如果create_option是conf::None，则返回false。这可能表示发生了错误条件。
// 如果是conf::Random，则调用creation::createRandom来创建与类型和空位相关的东西。
// 如果是conf::Pipe，则调用creation::createFromPile并记录一条消息。
// 如果是conf::Restart，则不执行任何操作。
// 最后，如果没有出现问题，函数返回true。
// 这段代码似乎是一个更大程序的一部分，该程序基于一些配置值设置了一个模拟环境，可能用于蒙特卡洛模拟。根据配置，它创建不同类型的晶格或初始化其他模拟参数。

bool PKMC::prepare() {
  // struct timeval time_start, time_end;
  // double timeuse;
  // gettimeofday(&time_start,NULL);
  mpi_types::setInterMPIType();
  sim = new simulation();
  // 获取config文件中的设置的参数
  conf::ConfigValues config_v = ConfigParsing::getInstance()->configValues;
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " getInstance ending !!!.\n");
  // if the seed is 0, reset the seed using random device.
  r::seedRelocate(&config_v.seeds.create_types, &config_v.seeds.create_vacancy, &config_v.seeds.event_selection,
                  &config_v.seeds.time_inc);
  sim->createDomain(config_v.box_size, config_v.lattice_const, config_v.cutoff_radius);
  sim->createLattice();
  // initialize lattices types. 初始化晶格类型。
  // for (auto it = config_v.create.types_ratio.begin(); it != config_v.create.types_ratio.end(); ++it) {
  //   std::cout << *it << " ";
  // }
  switch (config_v.create.create_option) {
  case conf::None:
    // todo log error
    return false;
  // 调用这里的
  case conf::Random:
    if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " createRandom starting !!!.\n");
    creation::createRandom(config_v.seeds.create_types, config_v.seeds.create_vacancy, sim->box->lattice_list,
                           sim->box->va_list, config_v.create.types, config_v.create.types_ratio,
                           config_v.create.va_count, sim->_p_domain);
    if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " createRandom ending !!!.\n");
    break;
  case conf::Pipe:
    creation::createFromPile(config_v.create.pipe_input_box, config_v.seeds.create_types, sim->box->lattice_list,
                             sim->box->va_list, config_v.create.types, config_v.create.types_ratio, sim->_p_domain);
    kiwi::logs::i("create", "placed {} defects from pipe file.\n", sim->box->va_list->mp.size());
    break;
  case conf::Restart:
    break;
  }
  //std::cout << " _lattices[21][71][4].id is : " << sim->box->lattice_list->_lattices[21][71][4].id << " _lattices[21][71][3].id is : " << sim->box->lattice_list->_lattices[21][71][3].id << std::endl;
  // gettimeofday(&time_end,NULL);
  // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0e3 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e3;
  // std::cout << " rank : " << kiwi::mpiUtils::global_process.own_rank << " prepare time : " << timeuse << " ms " << std::endl;
  return true;
}

void PKMC::onStart() {
  // struct timeval time_start, time_end;
  // double timeuse;
  // gettimeofday(&time_start,NULL);
  //std::cout << " _)))lattices[21][71][4].id is : " << sim->box->lattice_list->_lattices[21][71][4].id << " _)))lattices[21][71][3].id is : " << sim->box->lattice_list->_lattices[21][71][3].id << std::endl;
  conf::ConfigValues config_v = ConfigParsing::getInstance()->configValues;
  //  set up ghost.
  //std::cout << " 11111111111111111111111111111111 " << std::endl;
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " prepareForStart starting !!!.\n");
  sim->prepareForStart();
  // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " prepareForStart ending !!!.\n");
  //std::cout << " _)))lattices[21][71][4].id is : " << sim->box->lattice_list->_lattices[21][71][4].id << " _)))lattices[21][71][3].id is : " << sim->box->lattice_list->_lattices[21][71][3].id << std::endl;
  // setup model and run simulation
  // 统计每个进程父负责的区域内的原子类型及其数量
  counter m_counter = counter::newCounter(sim->box->lattice_list); // atoms counter
  m_counter.setLatTypeToStrFunc(&lat::LatTypesString);
  std::cout << m_counter;
  ABVIModel model(sim->box, config_v.attempt_freq, config_v.temperature);
  model.setColoredDomain(sim->_p_domain);
  sim->get_neighbour_local_sub_box();
  initialize_gpu(kiwi::mpiUtils::global_process.own_rank);
  // hipmalloc并且分配初始化一些常量
  // gpu_prepare(config_v.attempt_freq, config_v.temperature, sim->box->lattice_list, sim->_p_domain);
  //gpu_prepare(config_v.attempt_freq, config_v.temperature, sim->box->lattice_list, sim->_p_domain);
  gpu_prepare(config_v.attempt_freq, config_v.temperature);
  MEventListener m_listener(m_counter, sim->box->lattice_list->meta);
  model.setEventListener(&m_listener);
  // run simulation
  MEventHook m_event_hook(config_v.output, sim->box->lattice_list, &m_counter);
  
  sim->simulate(&model, &m_event_hook, config_v.seeds.time_inc, config_v.physics_time);
  // gettimeofday(&time_end,NULL);
  // timeuse = (time_end.tv_sec - time_start.tv_sec) * 1.0e3 + (double)(time_end.tv_usec - time_start.tv_usec) / 1.0e3;
  // std::cout << " rank : " << kiwi::mpiUtils::global_process.own_rank << " onStart time : " << timeuse << " ms " << std::endl;
}

void PKMC::onFinish() { mpi_types::unsetInterMPIType(); }

void PKMC::beforeDestroy() { delete sim; }

void PKMC::onDestroy() { kiwiApp::onDestroy(); }
