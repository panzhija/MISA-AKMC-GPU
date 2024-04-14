//
// Created by genshen on 2019/12/10.
//

#ifndef MISA_KMC_SEED_RELOCATE_H
#define MISA_KMC_SEED_RELOCATE_H

#include <cstdint>
#include <logs/logs.h>
#include <random>
#include <utils/random/random.h>

namespace r {
  /**
   * \brief if the seed is 0, set the seed using random device.
   */
// 这段代码定义了一个函数 seedRelocate，该函数用于重新定位随机数生成的种子值。它接受四个指向 uint32_t 类型的指针作为参数，分别用于四个不同的种子值：create_types、create_vacancy、event_selection 和 time_inc。
// 函数的主要作用是检查这些种子值是否等于特殊值 r::seed_auto，如果等于则使用 std::random_device 生成一个随机种子值，并将其赋值给相应的参数。
// 然后，函数会记录生成的种子值到日志中，以便追踪种子的值。
// 这个函数的用途通常是在随机数生成过程中，如果用户没有显式指定种子值，那么就使用真正的随机数来初始化这些种子值，以确保每次运行程序时都会得到不同的随机数序列。
// 这可以增加模拟程序的随机性和可重复性，因为用户可以通过设置相同的种子值来获得相同的模拟结果。
// 需要注意的是，std::random_device 通常用于生成真正的随机数，因此生成的种子值会受到系统随机源的影响，这可以提高随机性。
// 这个函数的实现假设 r::seed_auto 是一个特殊的标记值，用于表示需要生成随机种子值。如果 r::seed_auto 的定义不在这段代码中，那么它可能在程序的其他地方定义。
  void seedRelocate(int64_t *create_types, int64_t *create_vacancy, int64_t *event_selection, int64_t *time_inc) {
    std::random_device rd;
    if (*create_types == r::seed_auto) {
      *create_types = rd();
      kiwi::logs::i("random", "create_types seed: {}\n", *create_types);
    }
    if (*create_vacancy == r::seed_auto) {
      *create_vacancy = rd();
      kiwi::logs::i("random", "create_vacancy seed: {}\n", *create_vacancy);
    }
    if (*event_selection == r::seed_auto) {
      *event_selection = rd();
      kiwi::logs::i("random", "event_selection seed: {}\n", *event_selection);
    }
    if (*time_inc == r::seed_auto) {
      *time_inc = rd();
      kiwi::logs::i("random", "time_inc seed: {}\n", *time_inc);
    }
  }
} // namespace r

#endif // MISA_KMC_SEED_RELOCATE_H
