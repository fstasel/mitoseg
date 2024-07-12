/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "phase1_main.h"

sem_t sem1;
pthread_mutex_t roilock = PTHREAD_MUTEX_INITIALIZER;

void mainFunctionPhase1(int s) {
    cv::Mat originalImage, roiImage;
    cv::Mat grayImage, alImage, resmpImage, visualResmpImage,
        visualResmpImage_color;
    cv::Mat bilateral, smoothed, visualBilateral, visualSmoothed;
    vector<ridge> r;
    cv::Mat ridges, visualRidges, ridgeMask;
    cv::Mat visualRidges_dir;
    vector<localEnergy> e_lo, e_hi;
    cv::Mat visualLE_lo, visualLE_hi;
    int lew_lo, leh_lo, lew_hi, leh_hi;
    vector<curve> clist_lo, clist_hi;
    vector<curve> clist_lo_flt, clist_hi_flt, clist_hix;
    cv::Mat visualCurves_lo, visualCurves_hi;
    cv::Mat visualCurves_hilo;

    // Load Image Slice
    originalImage = loadSlice(s);

    // Get image roi
    if (s != SLICE_START)
        pthread_mutex_lock(&roilock);
    // roi image is created inside function as type of original image
    roiImage = getImageROI(originalImage);
    pthread_mutex_unlock(&roilock);

    // Convert to gray image
    cv::cvtColor(roiImage, grayImage, cv::COLOR_BGR2GRAY);

    // Auto levels
    alImage.create(roiImage.size(), CV_8UC1);
    autoLevels(grayImage, alImage);

    // Resample image; created inside function as CV_32FC1
    resmpImage = getResampledImage(alImage);

    // Create Core Images
    bilateral.create(resmpImage.size(), CV_32FC1);
    smoothed.create(resmpImage.size(), CV_32FC1);
    r.resize(resmpImage.cols * resmpImage.rows);
    ridges.create(resmpImage.size(), CV_32FC1);
    ridgeMask.create(resmpImage.size(), CV_8UC1);

    // Create Visual Images
    visualResmpImage.create(resmpImage.size(), CV_8UC1);
    visualResmpImage_color.create(resmpImage.size(), CV_8UC3);
    visualBilateral.create(resmpImage.size(), CV_8UC1);
    visualSmoothed.create(resmpImage.size(), CV_8UC1);
    visualRidges.create(resmpImage.size(), CV_8UC1);
    visualRidges_dir.create(resmpImage.size(), CV_8UC3);
    // visualLE_lo = created inside function...
    // visualLE_hi = created inside function...
    visualCurves_lo.create(resmpImage.size(), CV_8UC3);
    visualCurves_hi.create(resmpImage.size(), CV_8UC3);
    visualCurves_hilo.create(resmpImage.size(), CV_8UC3);

    // Pre-processing
    // Smooth
    fastBilateralFilter(resmpImage, bilateral, (int)BSMOOTH_SPAT, BSMOOTH_GRAY,
                        0);

    cv::GaussianBlur(bilateral, smoothed, cv::Size(0, 0), SMOOTH_STDEV, 0.0,
                     cv::BORDER_REPLICATE);
    // ByPass:cvSmooth(resmpImage, smoothed, CV_GAUSSIAN, 0, 0, SMOOTH_STDEV, 0
    // );	//by-pass bilateral filtering
    /////////////////

    // Get Ridges and normalize (autolevels algorithm)
    getRidges(smoothed, r, ridges, ridgeMask);
    normalizeRidges(r, ridges, ridgeMask);

    // Get local energy
    e_lo = getLocalEnergyMap(r, ridges.cols, ridges.rows, (int)LE_SSIZE_LO,
                             (int)LE_BSIZE_LO, &lew_lo, &leh_lo);
    e_hi = getLocalEnergyMap(r, ridges.cols, ridges.rows, (int)LE_SSIZE_HI,
                             (int)LE_BSIZE_HI, &lew_hi, &leh_hi);

    // Fit curves
    initCurveFitting();
    clist_lo = fitCurves(e_lo, lew_lo, leh_lo, FCS_XYSTEP_LO, FCS_XYRANGE_LO,
                         FCS_INIT_THRESH_LO, FCS_MINLEN_LO);
    clist_hi = fitCurves(e_hi, lew_hi, leh_hi, FCS_XYSTEP_HI, FCS_XYRANGE_HI,
                         FCS_INIT_THRESH_HI, FCS_MINLEN_HI);

    // Filter curves
    clist_lo_flt =
        filterCurves(clist_lo, (float)FLC_THRESH_LO, (float)FLC_AVGTHRESH_LO);
    clist_hi_flt =
        filterCurves(clist_hi, (float)FLC_THRESH_HI, (float)FLC_AVGTHRESH_HI);
    clist_hix = getHixCurves(clist_lo_flt, lew_lo, leh_lo, clist_hi_flt);

    // Visualize
    cv::normalize(resmpImage, visualResmpImage, 0, 255, cv::NORM_MINMAX);
    std::vector<cv::Mat> channels(3, visualResmpImage);
    cv::merge(channels, visualResmpImage_color);
    cv::normalize(bilateral, visualBilateral, 255, 0, cv::NORM_MINMAX);
    cv::normalize(smoothed, visualSmoothed, 255, 0, cv::NORM_MINMAX);
    cv::normalize(ridges, visualRidges, 255, 0, cv::NORM_MINMAX);
    visualizeRidgesDir(4, r, visualRidges_dir);
    visualLE_lo = visualizeLEMajority(e_lo, lew_lo, leh_lo, VISUAL_LE_SSIZE, 1);
    visualLE_hi = visualizeLEMajority(e_hi, lew_hi, leh_hi, VISUAL_LE_SSIZE, 1);
    visualizeCurves(visualResmpImage, clist_lo, (int)LE_SSIZE_LO,
                    visualCurves_lo, (float)FLC_THRESH_LO,
                    (float)FLC_AVGTHRESH_LO);
    visualizeCurves(visualResmpImage, clist_hi, (int)LE_SSIZE_HI,
                    visualCurves_hi, (float)FLC_THRESH_HI,
                    (float)FLC_AVGTHRESH_HI);
    visualizeHiLoCurves(clist_lo_flt, lew_lo, leh_lo, clist_hi_flt,
                        CURV_HILO_COV_TH, visualResmpImage, visualCurves_hilo);

    // Save data
    saveImageDat(s, "resmpImage_", resmpImage);
    saveRidgeStructure(s, "ridgeStruc_", r);
    saveLocalEnergyStructure(s, "e_lo_", e_lo, lew_lo, leh_lo);
    saveLocalEnergyStructure(s, "e_hi_", e_hi, lew_hi, leh_hi);
    saveCurveStructure(s, "clist_lo_", clist_lo);
    saveCurveStructure(s, "clist_hi_", clist_hi);
    saveCurveStructure(s, "clist_lo_flt_", clist_lo_flt);
    saveCurveStructure(s, "clist_hi_flt_", clist_hi_flt);
    saveCurveStructure(s, "clist_hix_", clist_hix);

    // Save Visual Images
    saveImage(s, "00_image_", "", visualResmpImage);
    saveImage(s, "01_bilateral_", "", visualBilateral);
    saveImage(s, "02_smoothed_", "", visualSmoothed);
    saveImage(s, "03_ridges_", "", visualRidges);
    saveImage(s, "04_ridges_dir_", "", visualRidges_dir);
    saveImage(s, "05_le_lo_", "", visualLE_lo);
    saveImage(s, "06_le_hi_", "", visualLE_hi);
    saveImage(s, "07_curves_lo_", "", visualCurves_lo);
    saveImage(s, "08_curves_hi_", "", visualCurves_hi);
    saveImage(s, "09_curves_hilo_", "", visualCurves_hilo);
}

void *threadPhase1(void *t_data) {
    int s = *(int *)t_data;

    cout << "Processing slice #" << s << "..." << endl;
    mainFunctionPhase1(s);

    sem_post(&sem1);
    pthread_exit(nullptr);
}

void startPhase1() {
    int s, t;
    int rc;
    int numThreads = SLICE_END - SLICE_START + 1;
    pthread_t *threads = new pthread_t[numThreads];
    pthread_attr_t attr;
    int *t_data = new int[numThreads];
    sem_init(&sem1, 0, numCores);
    pthread_mutex_trylock(&roilock);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    cout << ">>>> PHASE #1" << endl;
    for (s = SLICE_START, t = 0; s <= SLICE_END; s++, t++) {
        sem_wait(&sem1);
        t_data[t] = s;
        rc = pthread_create(&threads[t], &attr, threadPhase1,
                            (void *)&t_data[t]);
        if (rc) {
            cerr << "Error: Unable to create thread!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    pthread_attr_destroy(&attr);
    for (t = 0; t < numThreads; t++) {
        rc = pthread_join(threads[t], nullptr);
        if (rc) {
            cerr << "Error: Unable to join thread!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    delete[] threads;
    delete[] t_data;
}
