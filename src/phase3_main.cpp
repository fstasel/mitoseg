/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "phase3_main.h"

void revalidateSnakes(vector<snake25d> &snakeArray, int n) {
    // Visualizer
    Visualizer v;
    v.setStartSlice(SLICE_START);
    v.setNumSlices(sn25d_t);
    v.setFileNameTag(DESTPATH, "00_image_", "");
    v.update();

    // w, h
    int w_lo, h_lo;
    int w_hi, h_hi;
    getMapSize(v.getFrameSize(), (int)LE_SSIZE_LO, &w_lo, &h_lo);
    getMapSize(v.getFrameSize(), (int)LE_SSIZE_HI, &w_hi, &h_hi);

    // Load & process curve data
    dataPacket dp;
    dp.cs_lo = loadCurveStack(SLICE_START, SLICE_END, "clist_lo_flt_");
    dp.cs_hi = loadCurveStack(SLICE_START, SLICE_END, "clist_hix_");
    dp.ces_lo = createCurveEnergyStack(dp.cs_lo, w_lo, h_lo);
    dp.ces_hi = createCurveEnergyStack(dp.cs_hi, w_hi, h_hi);
    dp.grad_ces_lo = createCurveEnergyGradientStack(dp.ces_lo);
    dp.cps_lo = createCurvePointsStack(dp.cs_lo);
    dp.cps_hi = createCurvePointsStack(dp.cs_hi, w_hi, h_hi);

    for (int i = 0; i < n; i++) {
        validateSnake25d(dp, snakeArray[i]);
        cout << "\rRevalidating... " << i + 1 << " / "
             << " ";
        flush(cout);
    }
    cout << endl;
}

void startPhase3() {
    int start_z = 0, end_z = 0, s;
    Visualizer v;

    cout << ">>>> PHASE #3" << endl;

    // Visualize initial points
    vector<cv::Point> pts;
    int t_n_pts = 0;
    for (s = SLICE_START; s + sn25d_t - 1 <= SLICE_END; s += sn25d_t) {
        start_z = s;
        end_z = s + sn25d_t - 1;
        // Load initial pts
        pts = loadPoints(start_z, end_z);
        t_n_pts += (sn25d_t * pts.size());
        // Visualizer
        v.setStartSlice(start_z);
        v.setNumSlices(1);
        v.setFileNameTag(DESTPATH, "00_image_", "");
        v.setMarkerBuf(pts, start_z, LE_SSIZE_LO);
        v.update();
        saveImage(start_z, "10_initPts_", "", v.outputImage);
    }
    v.clearMarkerBuf();

    // Load snakes
    int z1 = SLICE_START;
    int z2 = end_z;
    vector<snake25d> snakeArray = loadSnakeArray(z1, z2);
    cout << snakeArray.size() << " snakes loaded..." << endl;
    // revalidateSnakes(snakeArray, n_snk);

    // Get valid polygons
    vector<poly25d> polyArray =
        convertValidSnakesToPolyArray(snakeArray, POLY_VALIDITY - EPS);
    cout << polyArray.size() << " valid polygons extracted..." << endl;

    // Get merged polygons
    vector<poly25d> polyArrayMerged =
        mergeArrayOfPoly25d(polyArray, POLY_MERGE);
    cout << polyArrayMerged.size() << " merged polygons extracted..." << endl;

    // Visualize polygons
    v.setFileNameTag(DESTPATH, "00_image_", "");
    v.setPolyBuf(polyArray, LE_SSIZE_LO);
    for (s = z1; s <= z2; s++) {
        v.setStartSlice(s);
        v.setNumSlices(1);
        v.update();
        saveImage(s, "11_mitos_", "", v.outputImage);
    }
    v.setPolyBuf(polyArrayMerged, LE_SSIZE_LO);
    for (s = z1; s <= z2; s++) {
        v.setStartSlice(s);
        v.setNumSlices(1);
        v.update();
        saveImage(s, "12_mitos_merged_", "", v.outputImage);
    }
    savePolyArrayAsPLY(polyArrayMerged);

    // Save as IMOD model
    cv::Mat img = loadSlice(SLICE_START);
    savePolyArrayAsIMOD(polyArrayMerged, img.cols, img.rows, SLICE_END,
                        LE_SSIZE_LO * TFACTOR / RESOLUTION, 1.0, ROI_X, ROI_Y);
}
