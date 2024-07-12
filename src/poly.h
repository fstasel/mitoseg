/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef POLY_H_
#define POLY_H_

#include <boost/format.hpp>
#include <opencv2/opencv.hpp>

#include "snake25d.h"

typedef vector<cv::Point> polygon;
typedef vector<polygon> poly2d;
struct poly25d {
    vector<poly2d> slice;
    uint t;
    int start_z, end_z;
};
typedef struct {
    int index[3];
} int3;

void initPoly2d(double (*)[2], int, poly2d &);
void initBlankPoly25d(int, int, poly25d &);
void deinitPoly25d(poly25d &);
void initPoly25d(snake25d &, poly25d &);
vector<poly25d> initPoly25dArray(snake25dList &);
void deinitPoly25dArray(vector<poly25d> &);
int poly2dArea(polygon &, int * = nullptr);
void correctPoly2d(poly2d &);
void combinePoly2d(const poly2d &, const poly2d &, int, poly2d &,
                   int * = nullptr, int * = nullptr, int * = nullptr);
void combinePoly25d(const poly25d &, const poly25d &, int, poly25d &,
                    int * = nullptr, int * = nullptr, int * = nullptr);
vector<poly25d> mergeArrayOfPoly25d(vector<poly25d> &, double);
vector<poly25d> convertValidSnakesToPolyArray(snake25dList &, double);
int getTriangles(polygon &, polygon &, vector<int3> &);
void savePolyArrayAsPLY(vector<poly25d> &);

#endif
