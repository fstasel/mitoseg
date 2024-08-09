/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CURVE_H_
#define CURVE_H_

#include <opencv2/opencv.hpp>

#include "preprocess.h"
#include "settings_parser.h"
using namespace mitoseg_settings;

struct localEnergy {
    vector<float> energy;
    int major;
    localEnergy() : energy(LE_BINS) {}
};

struct curve {
    // Curve fitness params
    float score, avgscore;
    int len;

    // Geometric model params
    int x1, y1, x2, y2, h;

    // Algebraic model params
    double R, a, b, sint, cost;
};

typedef vector<curve> curveList;

template <typename T> struct coord {
    T x;
    T y;
};
template <typename T> using curvePoints = vector<coord<T>>;

struct curveStack {
    vector<curveList> cstack;
    int start_z, end_z, t;
};

template <typename T> struct curvePointsStack {
    vector<curvePoints<T>> cps;
    int start_z, end_z, t;
};

vector<localEnergy> getLocalEnergyMap(vector<ridge> &, const int, const int,
                                      const int, const int, int *, int *);
cv::Mat visualizeLEMajority(vector<localEnergy> &, const int, const int,
                            const int, const int);
void calcCurveParams(vector<localEnergy> &, const int, const int, curve &,
                     const double);
void drawBinaryCurve(cv::Mat &, curve &, const int, const int, const double,
                     cv::Scalar = cv::Scalar(255, 0, 0));
void drawCurve(cv::Mat &, curve &, const int, const double, cv::Scalar);
int fitCurve(vector<localEnergy> &, const int, const int, curve &, const int,
             char * = nullptr);
float energystr(localEnergy &);
curveList fitCurves(vector<localEnergy> &, const int, const int, const int,
                    const int, const float, const int);
void visualizeCurves(cv::Mat &, curveList &, const int, cv::Mat &,
                     const float = -1.0f, const float = -1.0f);
void initCurveFitting();
cv::Mat drawBinaryCurves(curveList &, const int, const int, const int);
void visualizeHiLoCurves(curveList &, const int, const int, curveList &,
                         const double, cv::Mat &, cv::Mat &);
curveList getHixCurves(curveList &, const int, const int, curveList &);
curveList filterCurves(const curveList &, const float = -1.0f,
                       const float = -1.0f);
curvePoints<double> getCurvePoints(const curve &, const double);
void saveLocalEnergyStructure(int, const string, vector<localEnergy> &, int,
                              int);
void loadLocalEnergyStructure(int, const string, vector<localEnergy> &, int *,
                              int *);
void saveCurveStructure(int, const string, curveList &);
curveList loadCurveStructure(int, const string);
int isInside(curve &, double, double, double = 0.0);

#endif
