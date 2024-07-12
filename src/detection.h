/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef DETECTION_H_
#define DETECTION_H_

#include <boost/format.hpp>
#include <opencv2/opencv.hpp>

#include "snake25d.h"

template <typename T> T quick_select(T[], int);
int isValidMitoSlice_25d(snake25d &, int, dataPacket &);

#endif
