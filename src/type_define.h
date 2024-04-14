//
// Created by genshen on 2018/11/9.
//

#ifndef MISA_KMC_TYPE_DEFINE_H
#define MISA_KMC_TYPE_DEFINE_H

#include "building_config.h"
#include <stdint.h>

/**
 * the type of lattice position in cartesian coordinate system
 */
typedef int64_t _type_box_size;
typedef int64_t _type_lattice_size;
typedef int64_t _type_lattice_offset;
typedef int64_t _type_lattice_coord;
typedef int64_t _type_lattice_id;
typedef int64_t _type_lattice_count;

/*!
 * \brief the type of transition rate
 */
typedef double _type_rate;

/**
 * \brief neighbor types
 */
typedef unsigned char _type_dirs_status; // todo compatible with _type_neighbour_status
typedef unsigned char _type_dir_id;
typedef unsigned char _type_dirs_size;

#endif // MISA_KMC_TYPE_DEFINE_H
