/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef PHASE1_MAIN_H_
#define PHASE1_MAIN_H_

#include <pthread.h>
#include <semaphore.h>

#include <opencv2/opencv.hpp>

#include "curve.h"
#include "settings_parser.h"
using namespace mitoseg_settings;

void mainFunctionPhase1(int);
void *threadPhase1(void *);
void startPhase1();

#endif
