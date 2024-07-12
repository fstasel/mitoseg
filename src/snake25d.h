/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SNAKE25D_H_
#define SNAKE25D_H_

#include <pthread.h>
#include <semaphore.h>

#include "curve.h"
#include "settings_parser.h"
using namespace mitoseg_settings;

enum SNAKE_STATE {
    SNAKE_GROWING,
    SNAKE_CONVERGED,
    SNAKE_MAXITER,
    SNAKE_MAXAREA,
    SNAKE_MINAREA
};

#define SNAKE_INVALID false
#define SNAKE_VALID true

struct snake25d {
    double node[SN25D_T][SN_N][2];
    double d1_xy[SN25D_T][SN_N][2];
    double d2_xy[SN25D_T][SN_N][2];
    double d4_xy[SN25D_T][SN_N][2];
    double d1_z[SN25D_T][SN_N][2];
    double d2_z[SN25D_T][SN_N][2];
    double d4_z[SN25D_T][SN_N][2];
    double grad_Eint[SN25D_T][SN_N][2];
    double grad_Eext[SN25D_T][SN_N][2];
    double area[SN25D_T];
    SNAKE_STATE status;
    bool isValid[SN25D_T];
    double validity;
    int start_z, end_z, t;
    double aspectRatio_z;
};

typedef vector<snake25d> snake25dList;

struct curveEnergyMap {
    vector<double> ce;
    int w, h;
};

struct curveEnergyStack {
    vector<curveEnergyMap> ces;
    int start_z, end_z, t;
};

struct gradVec {
    double x, y;
};

struct curveEnergyGradientMap {
    vector<gradVec> grad_ce;
    int w, h;
};

struct curveEnergyGradientStack {
    vector<curveEnergyGradientMap> grad_ces;
    int start_z, end_z, t;
};

struct dataPacket {
    curveStack cs_lo;
    curveStack cs_hi;
    curveEnergyStack ces_lo;
    curveEnergyStack ces_hi;
    curveEnergyGradientStack grad_ces_lo;
    curvePointsStack<double> cps_lo;
    curvePointsStack<int> cps_hi;
};

struct t_param {
    dataPacket *dp;
    double ix, iy;
    int start_z, end_z;
    double ar_z;
    snake25dList *outputSnakeList;
    volatile int *returnValue;
};
//
curveStack loadCurveStack(int, int, const string &);
curvePoints<double> createCurvePoints(curveList &);
curvePoints<int> createCurvePoints(curveList &, int, int);
curvePointsStack<double> createCurvePointsStack(curveStack &);
curvePointsStack<int> createCurvePointsStack(curveStack &, int, int);
curveEnergyMap createCurveEnergyMap(curveList &, const int, const int);
curveEnergyStack createCurveEnergyStack(curveStack &, const int, const int);
curveEnergyGradientMap createCurveEnergyGradientMap(curveEnergyMap &);
curveEnergyGradientStack createCurveEnergyGradientStack(curveEnergyStack &);
//
void getMapSize(cv::Size, const int, int *, int *);
void compSnakeSliceArea(snake25d &, const int);
void compSnakeArea(snake25d &);
void compApproxIntersectionArea(double[][2], int, double[][2], int, double *,
                                double *, double *, double = 1.0);
double getDiceSimilarity2d(double[][2], int, double[][2], int, double = 1.0);
double getDiceSimilaritySnake25d(snake25d &, snake25d &, double = 1.0);
snake25d initSnake25d(const double, const double, const int, const int,
                      const double, const int, const int, const double);
void compSnakeDerivatives(snake25d &);
void compIntEnergyGrad(snake25d &, const double, const double, const double,
                       const double);
void compExtEnergyGrad(snake25d &, curveEnergyGradientStack &, const double,
                       const double);
void equidistantCorrectionUV(snake25d &, double[][SN_N][2], double[][SN_N][2]);
void compUpdateVectors(snake25d &, double[][SN_N][2]);
double alphaCoefficient(double);
void updateSnake(snake25d &, double[][SN_N][2], const int, const int,
                 double[][SN_N][2], double *);
void resetMovement(double[][SN_N][2], double *);
snake25d startSnake25d(curveEnergyGradientStack &, const double, const double,
                       const double, const double, const double, const double,
                       const double, snake25d *, const double, const double,
                       const int, int);
int inflationSnake25d(dataPacket &, const double, const double, const int, int,
                      const double, const double, const double, const double,
                      const double, const double, const double, const double,
                      const double, snake25dList &);
void reportValidation(dataPacket &, snake25d &);
void validateSnake25d(dataPacket &, snake25d &);
vector<cv::Point> findInitialPointsWithinZRange(dataPacket &, int, int);
void *threadPhase2(void *);
void retrieveAllSnakes25d(dataPacket &, vector<cv::Point> &, int, int,
                          snake25dList &);
void saveSnakeList(snake25dList &, int, int);
snake25dList loadSnakeArray(int, int);
void savePoints(vector<cv::Point> &, int, int);
vector<cv::Point> loadPoints(int, int);
snake25dList filterSnakeListByValidity(snake25dList &, double);
vector<gradVec> getEnergyGradient(const vector<double> &, const int, const int);
vector<double> getCurveEnergyImage(vector<curve> &, const int, const int);

#endif
