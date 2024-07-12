/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "phase2_main.h"

void startPhase2() {
    cout << ">>>> PHASE #2" << endl;

    // Initial Points
    vector<cv::Point> pts;

    // Snakes
    snake25dList snakeList;

    // w, h
    int w_lo, h_lo;
    int w_hi, h_hi;
    cv::Mat rsmpImg = loadImageDat(SLICE_START, "resmpImage_");
    getMapSize(cv::Size(rsmpImg.cols, rsmpImg.rows), (int)LE_SSIZE_LO, &w_lo,
               &h_lo);
    getMapSize(cv::Size(rsmpImg.cols, rsmpImg.rows), (int)LE_SSIZE_HI, &w_hi,
               &h_hi);

    // Load & process curve data
    dataPacket dp;
    dp.cs_lo = loadCurveStack(SLICE_START, SLICE_END, "clist_lo_flt_");
    dp.cs_hi = loadCurveStack(SLICE_START, SLICE_END, "clist_hix_");
    dp.ces_lo = createCurveEnergyStack(dp.cs_lo, w_lo, h_lo);
    dp.ces_hi = createCurveEnergyStack(dp.cs_hi, w_hi, h_hi);
    dp.grad_ces_lo = createCurveEnergyGradientStack(dp.ces_lo);
    dp.cps_lo = createCurvePointsStack(dp.cs_lo);
    dp.cps_hi = createCurvePointsStack(dp.cs_hi, w_hi, h_hi);

    int s, start_z, end_z = 0;

    for (s = SLICE_START; s + sn25d_t - 1 <= SLICE_END; s += sn25d_t) {
        start_z = s;
        end_z = s + sn25d_t - 1;

        pts = findInitialPointsWithinZRange(dp, start_z, end_z);

        // Save data
        cout << "Saving " << pts.size() << " initial points..." << endl;
        savePoints(pts, start_z, end_z);

        retrieveAllSnakes25d(dp, pts, start_z, end_z, snakeList);
    }

    // Save data
    cout << "Saving " << snakeList.size() << " snakes..." << endl;
    saveSnakeList(snakeList, SLICE_START, end_z);
}
