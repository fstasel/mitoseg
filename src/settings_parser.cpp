/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "settings_parser.h"

#include <yaml-cpp/yaml.h>

#include <iostream>
using namespace std;

namespace mitoseg_settings {
    // General settings
    int numCores = 1;
    int THREAD_STACK_SIZE = 32;
    string SETTINGS_PATH = "";

    // Phase 1 - Visualization parameters
    int VISUAL_LE_SSIZE = 8;

    // Phase 1 - Preprocessing parameters
    double TFACTOR = 2.0;
    double AUTOLEVELSCUT = 0.005;
    double BSMOOTH_SPAT = 60 / TFACTOR;
    double BSMOOTH_GRAY = 0.2;
    double SMOOTH_STDEV = 3.0 / TFACTOR;

    // Phase 1 - Common curve fit parameters
    double FC_MAXCURV = 1.0;
    int FC_MAXITER = 200;
    double FC_STEP = 1.0;
    int FC_XYRANGE = 2;
    int FC_XYSTEP = 1;
    int FC_HRANGE = 2;
    int FC_HSTEP = 1;
    double FC_TOL = 0.01;
    int LE_BINS = 4;

    // Phase 1 - Low frequency curve fit parameters
    double LE_SSIZE_LO = 4 / TFACTOR;
    double LE_BSIZE_LO = 30 / TFACTOR;
    double FCS_XYSTEP_LO = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
    double FCS_XYRANGE_LO = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
    double FCS_MINLEN_LO = ((int)((100 / TFACTOR) / LE_SSIZE_LO));
    double FLC_THRESH_LO = (100.0 * LE_BSIZE_LO / LE_SSIZE_LO);
    double FLC_AVGTHRESH_LO = (1.2 * LE_BSIZE_LO);
    double FCS_INIT_THRESH_LO = 0.4;

    // Phase 1 - High frequency curve fit parameters
    double LE_SSIZE_HI = (4 / TFACTOR);
    double LE_BSIZE_HI = (8 / TFACTOR);
    double FCS_XYSTEP_HI = ((int)((16 / TFACTOR) / LE_SSIZE_HI));
    double FCS_XYRANGE_HI = ((int)((16 / TFACTOR) / LE_SSIZE_HI));
    double FCS_MINLEN_HI = ((int)((20 / TFACTOR) / LE_SSIZE_HI));
    double FLC_THRESH_HI = (20 * LE_BSIZE_HI / LE_SSIZE_HI);
    double FLC_AVGTHRESH_HI = (0.7 * LE_BSIZE_HI);
    double FCS_INIT_THRESH_HI = 0.4;

    // Phase 1 - Low frequency curves coverage distance
    double CURV_HILO_COV = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
    double CURV_HILO_COV_TH = 0.7;

    // Phase 1 - Detection parameters
    double DT_BOUNDARY_T = 20.0;
    double DT_REGION_T = 0.1;
    double DT_GAP_T = 15.0;
    double DT_GAPTOTAL_TD = 600.0;
    double DT_GAPTOTAL_T = ((DT_GAPTOTAL_TD / TFACTOR) / LE_SSIZE_LO);
    double DT_GAPMAX_TD = 600.0;
    double DT_GAPMAX_T = ((DT_GAPMAX_TD / TFACTOR) / LE_SSIZE_LO);
    double DT_GAPTOTALRATIO_T = 0.4;
    double DT_GAPBORDERRATIO_T = 0.4;
    int DT_GAPBORDERCOUNT_T = 1;

    double DT_CURVG_T = 15.0;
    double DT_CURVL_T = 0.6;
    double DT_CURVL_SEG = 0.25;
    double DT_SIGN_SMOOTH = 0.05;
    double DT_SIGN_RATIO = 4.0;
    int DT_SIGN_MAXNUMCRIT_T = 4;
    double DT_MAJORAXIS_T = ((2000.0 / TFACTOR) / LE_SSIZE_LO);
    double DT_MINORAXIS_T = ((140.0 / TFACTOR) / LE_SSIZE_LO);
    double DT_MINAXIS_T = ((70.0 / TFACTOR) / LE_SSIZE_LO);
    double DT_MINMINOR_RATIO_T = 0.2;
    double DT_MINAREA =
        ((20000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
    double DT_MAXAREA =
        ((700000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
    bool DT_REPORT = false;
    int DT25D_MED_RANGE = 5;

    // Phase 2 - Common snake parameters
    // int SN_N = 100;
    double SN_GAUSSIAN = 1.0;
    double SN_MAXAREA = DT_MAXAREA;
    double SN_MINAREA = DT_MINAREA;

    // Phase 2 - 2.5D snake parameters
    double SN25D_W_TENSION = 1.0;
    double SN25D_W_CURVATURE = 200.0;
    double SN25D_W_ZTENSION = 5.0;
    double SN25D_W_ZCURVATURE = 5.0;
    double SN25D_W_ECURVE = 1.0;
    double SN25D_W_EINF_MIN = 0.5;
    double SN25D_W_EINF_MAX = 3.0;
    double SN25D_W_EINF_STEP = 0.5;
    // int SN25D_T = 500;
    int sn25d_t = 20;
    double SN25D_K = 1.0;
    double SN25D_INITR = 10.0;
    int SN25D_MAXITER = 400000;
    int SN25D_SHORTTERM_CONV_ITER = 10000;
    double SN25D_SHORTTERM_CONV = 2.0;
    int SN25D_LONGTERM_CONV_ITER = 40000;
    double SN25D_LONGTERM_CONV = 10.0;
    double SN25D_MAXAREA =
        ((700000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
    double SN25D_MINAREA =
        ((1000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
    double SN25D_INF_CONV = 0.95;
    double SN25D_INITPTS_EPS = ((100 / TFACTOR) / LE_SSIZE_LO);
    double SN25D_INITPTS_MIN = ((int)(1.5 * sn25d_t));
    double SN25D_ALPHA = 0.0;

    // Phase 3 parameters
    double POLY_VALIDITY = 0.75;
    double POLY_MERGE = 0.3;

    // Phase 3 - Visualizer parameters
    int MAXIMAGES = 150;
    int MAXBUFFER = 150;
    int WINDOW_WIDTH = 720;
    int WINDOW_HEIGHT = 720;

    bool loadSettings() {
        YAML::Node settings;
        try {
            settings = YAML::LoadFile(SETTINGS_PATH);
        } catch (const exception &e) {
            return false;
        }

        if (settings["numCores"])
            numCores = settings["numCores"].as<int>();
        if (settings["THREAD_STACK_SIZE"])
            THREAD_STACK_SIZE = settings["THREAD_STACK_SIZE"].as<int>();
        if (settings["VISUAL_LE_SSIZE"])
            VISUAL_LE_SSIZE = settings["VISUAL_LE_SSIZE"].as<int>();
        if (settings["TFACTOR"])
            TFACTOR = settings["TFACTOR"].as<double>();
        if (settings["AUTOLEVELSCUT"])
            AUTOLEVELSCUT = settings["AUTOLEVELSCUT"].as<double>();
        BSMOOTH_SPAT = 60 / TFACTOR;
        if (settings["BSMOOTH_GRAY"])
            BSMOOTH_GRAY = settings["BSMOOTH_GRAY"].as<double>();
        SMOOTH_STDEV = 3.0 / TFACTOR;
        if (settings["FC_MAXCURV"])
            FC_MAXCURV = settings["FC_MAXCURV"].as<double>();
        if (settings["FC_MAXITER"])
            FC_MAXITER = settings["FC_MAXITER"].as<int>();
        if (settings["FC_STEP"])
            FC_STEP = settings["FC_STEP"].as<double>();
        if (settings["FC_XYRANGE"])
            FC_XYRANGE = settings["FC_XYRANGE"].as<int>();
        if (settings["FC_XYSTEP"])
            FC_XYSTEP = settings["FC_XYSTEP"].as<int>();
        if (settings["FC_HRANGE"])
            FC_HRANGE = settings["FC_HRANGE"].as<int>();
        if (settings["FC_HSTEP"])
            FC_HSTEP = settings["FC_HSTEP"].as<int>();
        if (settings["FC_TOL"])
            FC_TOL = settings["FC_TOL"].as<double>();
        if (settings["LE_BINS"])
            LE_BINS = settings["LE_BINS"].as<int>();
        LE_SSIZE_LO = 4 / TFACTOR;
        LE_BSIZE_LO = 30 / TFACTOR;
        FCS_XYSTEP_LO = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
        FCS_XYRANGE_LO = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
        FCS_MINLEN_LO = ((int)((100 / TFACTOR) / LE_SSIZE_LO));
        FLC_THRESH_LO = (100.0 * LE_BSIZE_LO / LE_SSIZE_LO);
        FLC_AVGTHRESH_LO = (1.2 * LE_BSIZE_LO);
        if (settings["FCS_INIT_THRESH_LO"])
            FCS_INIT_THRESH_LO = settings["FCS_INIT_THRESH_LO"].as<double>();
        LE_SSIZE_HI = (4 / TFACTOR);
        LE_BSIZE_HI = (8 / TFACTOR);
        FCS_XYSTEP_HI = ((int)((16 / TFACTOR) / LE_SSIZE_HI));
        FCS_XYRANGE_HI = ((int)((16 / TFACTOR) / LE_SSIZE_HI));
        FCS_MINLEN_HI = ((int)((20 / TFACTOR) / LE_SSIZE_HI));
        FLC_THRESH_HI = (20 * LE_BSIZE_HI / LE_SSIZE_HI);
        FLC_AVGTHRESH_HI = (0.7 * LE_BSIZE_HI);
        if (settings["FCS_INIT_THRESH_HI"])
            FCS_INIT_THRESH_HI = settings["FCS_INIT_THRESH_HI"].as<double>();
        CURV_HILO_COV = ((int)((60 / TFACTOR) / LE_SSIZE_LO));
        if (settings["CURV_HILO_COV_TH"])
            CURV_HILO_COV_TH = settings["CURV_HILO_COV_TH"].as<double>();
        if (settings["DT_BOUNDARY_T"])
            DT_BOUNDARY_T = settings["DT_BOUNDARY_T"].as<double>();
        if (settings["DT_REGION_T"])
            DT_REGION_T = settings["DT_REGION_T"].as<double>();
        if (settings["DT_GAP_T"])
            DT_GAP_T = settings["DT_GAP_T"].as<double>();
        if (settings["DT_GAPTOTAL_TD"])
            DT_GAPTOTAL_TD = settings["DT_GAPTOTAL_TD"].as<double>();
        DT_GAPTOTAL_T = ((DT_GAPTOTAL_TD / TFACTOR) / LE_SSIZE_LO);
        if (settings["DT_GAPMAX_TD"])
            DT_GAPMAX_TD = settings["DT_GAPMAX_TD"].as<double>();
        DT_GAPMAX_T = ((DT_GAPMAX_TD / TFACTOR) / LE_SSIZE_LO);
        if (settings["DT_GAPTOTALRATIO_T"])
            DT_GAPTOTALRATIO_T = settings["DT_GAPTOTALRATIO_T"].as<double>();
        if (settings["DT_GAPBORDERRATIO_T"])
            DT_GAPBORDERRATIO_T = settings["DT_GAPBORDERRATIO_T"].as<double>();
        if (settings["DT_GAPBORDERCOUNT_T"])
            DT_GAPBORDERCOUNT_T = settings["DT_GAPBORDERCOUNT_T"].as<int>();
        if (settings["DT_CURVG_T"])
            DT_CURVG_T = settings["DT_CURVG_T"].as<double>();
        if (settings["DT_CURVL_T"])
            DT_CURVL_T = settings["DT_CURVL_T"].as<double>();
        if (settings["DT_CURVL_SEG"])
            DT_CURVL_SEG = settings["DT_CURVL_SEG"].as<double>();
        if (settings["DT_SIGN_SMOOTH"])
            DT_SIGN_SMOOTH = settings["DT_SIGN_SMOOTH"].as<double>();
        if (settings["DT_SIGN_RATIO"])
            DT_SIGN_RATIO = settings["DT_SIGN_RATIO"].as<double>();
        if (settings["DT_SIGN_MAXNUMCRIT_T"])
            DT_SIGN_MAXNUMCRIT_T = settings["DT_SIGN_MAXNUMCRIT_T"].as<int>();
        DT_MAJORAXIS_T = ((2000.0 / TFACTOR) / LE_SSIZE_LO);
        DT_MINORAXIS_T = ((140.0 / TFACTOR) / LE_SSIZE_LO);
        DT_MINAXIS_T = ((70.0 / TFACTOR) / LE_SSIZE_LO);
        if (settings["DT_MINMINOR_RATIO_T"])
            DT_MINMINOR_RATIO_T = settings["DT_MINMINOR_RATIO_T"].as<double>();
        DT_MINAREA =
            ((20000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
        DT_MAXAREA =
            ((700000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
        if (settings["DT_REPORT"])
            DT_REPORT = settings["DT_REPORT"].as<bool>();
        if (settings["DT25D_MED_RANGE"])
            DT25D_MED_RANGE = settings["DT25D_MED_RANGE"].as<int>();
        // if (settings["SN_N"])
        //     SN_N = settings["SN_N"].as<int>();
        if (settings["SN_GAUSSIAN"])
            SN_GAUSSIAN = settings["SN_GAUSSIAN"].as<double>();
        SN_MAXAREA = DT_MAXAREA;
        SN_MINAREA = DT_MINAREA;
        if (settings["SN25D_W_TENSION"])
            SN25D_W_TENSION = settings["SN25D_W_TENSION"].as<double>();
        if (settings["SN25D_W_CURVATURE"])
            SN25D_W_CURVATURE = settings["SN25D_W_CURVATURE"].as<double>();
        if (settings["SN25D_W_ZTENSION"])
            SN25D_W_ZTENSION = settings["SN25D_W_ZTENSION"].as<double>();
        if (settings["SN25D_W_ZCURVATURE"])
            SN25D_W_ZCURVATURE = settings["SN25D_W_ZCURVATURE"].as<double>();
        if (settings["SN25D_W_ECURVE"])
            SN25D_W_ECURVE = settings["SN25D_W_ECURVE"].as<double>();
        if (settings["SN25D_W_EINF_MIN"])
            SN25D_W_EINF_MIN = settings["SN25D_W_EINF_MIN"].as<double>();
        if (settings["SN25D_W_EINF_MAX"])
            SN25D_W_EINF_MAX = settings["SN25D_W_EINF_MAX"].as<double>();
        if (settings["SN25D_W_EINF_STEP"])
            SN25D_W_EINF_STEP = settings["SN25D_W_EINF_STEP"].as<double>();
        // if (settings["SN25D_T"])
        //     SN25D_T = settings["SN25D_T"].as<int>();
        if (settings["SN25D_K"])
            SN25D_K = settings["SN25D_K"].as<double>();
        if (settings["SN25D_INITR"])
            SN25D_INITR = settings["SN25D_INITR"].as<double>();
        if (settings["SN25D_MAXITER"])
            SN25D_MAXITER = settings["SN25D_MAXITER"].as<double>();
        if (settings["SN25D_SHORTTERM_CONV_ITER"])
            SN25D_SHORTTERM_CONV_ITER =
                settings["SN25D_SHORTTERM_CONV_ITER"].as<int>();
        if (settings["SN25D_SHORTTERM_CONV"])
            SN25D_SHORTTERM_CONV =
                settings["SN25D_SHORTTERM_CONV"].as<double>();
        if (settings["SN25D_LONGTERM_CONV_ITER"])
            SN25D_LONGTERM_CONV_ITER =
                settings["SN25D_LONGTERM_CONV_ITER"].as<int>();
        if (settings["SN25D_LONGTERM_CONV"])
            SN25D_LONGTERM_CONV = settings["SN25D_LONGTERM_CONV"].as<double>();
        SN25D_MAXAREA =
            ((700000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
        SN25D_MINAREA =
            ((1000 / (TFACTOR * TFACTOR)) / (LE_SSIZE_LO * LE_SSIZE_LO));
        if (settings["SN25D_INF_CONV"])
            SN25D_INF_CONV = settings["SN25D_INF_CONV"].as<double>();
        SN25D_INITPTS_EPS = ((100 / TFACTOR) / LE_SSIZE_LO);
        SN25D_INITPTS_MIN = ((int)(1.5 * sn25d_t));
        if (settings["POLY_VALIDITY"])
            POLY_VALIDITY = settings["POLY_VALIDITY"].as<double>();
        if (settings["POLY_MERGE"])
            POLY_MERGE = settings["POLY_MERGE"].as<double>();
        if (settings["MAXIMAGES"])
            MAXIMAGES = settings["MAXIMAGES"].as<int>();
        if (settings["MAXBUFFER"])
            MAXBUFFER = settings["MAXBUFFER"].as<int>();
        if (settings["WINDOW_WIDTH"])
            WINDOW_WIDTH = settings["WINDOW_WIDTH"].as<int>();
        if (settings["WINDOW_HEIGHT"])
            WINDOW_HEIGHT = settings["WINDOW_HEIGHT"].as<int>();

        if (settings["SN25D_ALPHA"])
            SN25D_ALPHA = settings["SN25D_ALPHA"].as<double>();

        return true;
    }
} // namespace mitoseg_settings
