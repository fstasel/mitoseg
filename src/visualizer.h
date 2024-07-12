/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef VISUALIZER_H_
#define VISUALIZER_H_

#include <opencv2/opencv.hpp>

#include "poly.h"
#include "settings_parser.h"
using namespace mitoseg_settings;

class Visualizer {
  public:
    cv::Mat outputImage;
    //
    Visualizer();
    ~Visualizer();
    void setSnakeBuf(vector<snake25d> &, double = 1.0);
    void setMarkerBuf(vector<cv::Point> &, int = -1, double = 1.0);
    void clearMarkerBuf();
    void setPolyBuf(vector<poly25d> &, double = 1.0);
    void setStartSlice(int);
    int getStartSlice();
    void setNumSlices(int);
    int getNumSlices();
    cv::Size getFrameSize();
    cv::Size getGridSize();
    cv::Size getOutputSize();
    void setFileNameTag(const string &, const string &, const string &);
    void update();
    void flushBuffer();
    void getSliceCoor(cv::Point, int *, cv::Point &);

  private:
    int startSlice;
    int numSlices;
    int w, h, gx, gy, gw, gh, width, height;
    string fname_path, fname_tag1, fname_tag2;
    cv::Size sliceSize;
    //
    vector<cv::Mat> buf;
    vector<int> bufIndex;
    int nextFree;
    //
    vector<snake25d> snakeBuf;
    double snake_scale;
    //
    vector<cv::Point> markerBuf;
    double marker_scale;
    int marker_z;
    //
    vector<poly25d> polyBuf;
    double poly_scale;
    //
    int getBuffer(int, cv::Mat &);
    void insertBuffer(int, cv::Mat);
    int loadSliceImage(int, cv::Mat &);
    cv::Mat getSlice(int);
    cv::Mat getFrame(int, cv::Size);
    void getNATemplate(cv::Mat &);
    void drawFrame(int, cv::Mat &);
    void drawSnake(cv::Mat &, snake25d &, int);
    void drawMarker(cv::Mat &);
    void drawPoly(cv::Mat &, poly25d &, int);
};

void visualizeSnake(cv::Mat &, const double[][2], const double,
                    cv::Scalar = cv::Scalar(0, 255, 255), int = 3,
                    cv::Scalar = cv::Scalar(255, 255, 255));

#endif
