/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef IMOD_H_
#define IMOD_H_

#include <boost/format.hpp>

#include "poly.h"

#define IMOD_REVERSE_ENDIAN 0

int imod_endian_reverse(int x);
float imod_endian_reverse(float x);

#define IMOD_VALUE(x) (IMOD_REVERSE_ENDIAN ? (x) : imod_endian_reverse(x))

#define IMOD_FILE_ID 0x494d4f44
#define IMOD_VERSION_ID 0x56312e32
#define IMOD_IEOF_ID 0x49454f46
#define IMOD_OBJT_ID 0x4f424a54
#define IMOD_CONT_ID 0x434f4e54
#define IMOD_MESH_ID 0x4d455348

struct imod_model {
    char name[128];
    int xmax, ymax, zmax;
    int objsize;
    unsigned int flags;
    int drawmode, mousemode;
    int blacklevel, whitelevel;
    float xoffset, yoffset, zoffset;
    float xscale, yscale, zscale;
    int object, contour, point;
    int res;
    int thresh;
    float pixsize;
    int units;
    int csum;
    float alpha, beta, gamma;
};

struct imod_objt {
    char name[64];
    char extra[64];
    int contsize;
    unsigned int flags;
    int axis;
    int drawmode;
    float red, green, blue;
    int pdrawsize;
    unsigned char symbol;
    unsigned char symsize;
    unsigned char linewidth2;
    unsigned char linewidth;
    unsigned char linesty;
    unsigned char symflags;
    unsigned char sympad;
    unsigned char trans;
    int meshsize;
    int surfsize;
};

struct imod_cont {
    int psize;
    unsigned int flags;
    int time;
    int surf;
    // float[psize][3] pt
};

struct imod_mesh {
    int vsize, lsize;
    unsigned int flag;
    short time;
    short surf;
    // float[vsize][3] vert
    // int[lsize] list
};

void savePolyArrayAsIMOD(vector<poly25d> &, int, int, int, double, double, int,
                         int);

#endif
