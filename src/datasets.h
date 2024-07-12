/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef DATASETS_H_
#define DATASETS_H_

#include <boost/format.hpp>
#include <fstream>
#include <opencv2/opencv.hpp>

using namespace std;

extern string FNAME;
extern string DESTPATH;
extern string SOURCEPATH;
extern int SLICE_START;
extern int SLICE_END;
extern int ROI_X;
extern int ROI_Y;
extern int ROI_W;
extern int ROI_H;
extern double RESOLUTION;

cv::Mat loadSlice(int s);
void saveImage(int, const string &, const string &, const cv::Mat &);
void saveImageDat(int, const string &, const cv::Mat &);
cv::Mat loadImageDat(int, const string &);
double operator~(const cv::Vec3b &);
double getDev(const cv::Mat &);
bool checkRow(const cv::Mat &, const int, const double &);
bool checkColumn(const cv::Mat &, const int, const double &);
void autoROI(const cv::Mat &);
cv::Mat getImageROI(cv::Mat &);

#endif
