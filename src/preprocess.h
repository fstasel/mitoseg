/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef PREPROCESS_H_
#define PREPROCESS_H_

#include <opencv2/opencv.hpp>

#include "datasets.h"
#include "definitions.h"
#include "settings_parser.h"
using namespace mitoseg_settings;

// Ridge Structure
typedef struct {
    float angle;
    float dir;
    int dir4;
    int dir8;
    float mag;
    float vec[2];
    float val[2];
    bool mask;
} ridge;
//////////

////////////
void autoLevels(cv::Mat &, cv::Mat &);
cv::Mat getResampledImage(const cv::Mat &);
void getRidges(cv::Mat &, vector<ridge> &, cv::Mat &, cv::Mat &);
void visualizeRidgesDir(int, vector<ridge> &, cv::Mat &);
void normalizeRidges(vector<ridge> &, cv::Mat &, cv::Mat &);
void fastBilateralFilter(cv::Mat &, cv::Mat &, int, float, int);
void saveRidgeStructure(int, const string, vector<ridge> &);
vector<ridge> loadRidgeStructure(int, const string &);
/////////////

#endif
