/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SETTINGS_PARSER_H_
#define SETTINGS_PARSER_H_

#include <string>
using namespace std;

#define SN_N 100
#define SN25D_T 500

namespace mitoseg_settings {
    // General settings
    extern int numCores;
    extern int THREAD_STACK_SIZE;
    extern string SETTINGS_PATH;

    // Phase 1 - Visualization parameters
    extern int VISUAL_LE_SSIZE;

    // Phase 1 - Preprocessing parameters
    extern double TFACTOR;
    extern double AUTOLEVELSCUT;
    extern double BSMOOTH_SPAT;
    extern double BSMOOTH_GRAY;
    extern double SMOOTH_STDEV;

    // Phase 1 - Common curve fit parameters
    extern double FC_MAXCURV;
    extern int FC_MAXITER;
    extern double FC_STEP;
    extern int FC_XYRANGE;
    extern int FC_XYSTEP;
    extern int FC_HRANGE;
    extern int FC_HSTEP;
    extern int LE_BINS;

    // Phase 1 - Low frequency curve fit parameters
    extern double LE_SSIZE_LO;
    extern double LE_BSIZE_LO;
    extern double FCS_XYSTEP_LO;
    extern double FCS_XYRANGE_LO;
    extern double FCS_MINLEN_LO;
    extern double FLC_THRESH_LO;
    extern double FLC_AVGTHRESH_LO;
    extern double FCS_INIT_THRESH_LO;

    // Phase 1 - High frequency curve fit parameters
    extern double LE_SSIZE_HI;
    extern double LE_BSIZE_HI;
    extern double FCS_XYSTEP_HI;
    extern double FCS_XYRANGE_HI;
    extern double FCS_MINLEN_HI;
    extern double FLC_THRESH_HI;
    extern double FLC_AVGTHRESH_HI;
    extern double FCS_INIT_THRESH_HI;

    // Phase 1 - Low frequency curves coverage distance
    extern double CURV_HILO_COV;
    extern double CURV_HILO_COV_TH;

    // Phase 1 - Detection parameters
    extern double DT_BOUNDARY_T;
    extern double DT_REGION_T;
    extern double DT_GAP_T;
    extern double DT_GAPTOTAL_TD;
    extern double DT_GAPTOTAL_T;
    extern double DT_GAPMAX_TD;
    extern double DT_GAPMAX_T;
    extern double DT_GAPTOTALRATIO_T;
    extern double DT_GAPBORDERRATIO_T;
    extern int DT_GAPBORDERCOUNT_T;

    extern double DT_CURVG_T;
    extern double DT_CURVL_T;
    extern double DT_CURVL_SEG;
    extern double DT_SIGN_SMOOTH;
    extern double DT_SIGN_RATIO;
    extern int DT_SIGN_MAXNUMCRIT_T;
    extern double DT_MAJORAXIS_T;
    extern double DT_MINORAXIS_T;
    extern double DT_MINAXIS_T;
    extern double DT_MINMINOR_RATIO_T;
    extern double DT_MINAREA;
    extern double DT_MAXAREA;
    extern bool DT_REPORT;
    extern int DT25D_MED_RANGE;

    // Phase 2 - Common snake parameters
    extern int sn_n;
    extern double SN_GAUSSIAN;
    extern double SN_MAXAREA;
    extern double SN_MINAREA;

    // Phase 2 - 2.5D snake parameters
    extern double SN25D_W_TENSION;
    extern double SN25D_W_CURVATURE;
    extern double SN25D_W_ZTENSION;
    extern double SN25D_W_ZCURVATURE;
    extern double SN25D_W_ECURVE;
    extern double SN25D_W_EINF_MIN;
    extern double SN25D_W_EINF_MAX;
    extern double SN25D_W_EINF_STEP;
    extern int sn25d_t;
    extern double SN25D_K;
    extern double SN25D_INITR;
    extern int SN25D_MAXITER;
    extern int SN25D_SHORTTERM_CONV_ITER;
    extern double SN25D_SHORTTERM_CONV;
    extern int SN25D_LONGTERM_CONV_ITER;
    extern double SN25D_LONGTERM_CONV;
    extern double SN25D_MAXAREA;
    extern double SN25D_MINAREA;
    extern double SN25D_INF_CONV;
    extern double SN25D_INITPTS_EPS;
    extern double SN25D_INITPTS_MIN;

    // Phase 3
    extern double POLY_VALIDITY;
    extern double POLY_MERGE;

    // Phase 3 - Visualizer parameters
    extern int MAXIMAGES;
    extern int MAXBUFFER;
    extern int WINDOW_WIDTH;
    extern int WINDOW_HEIGHT;

    bool loadSettings();
} // namespace mitoseg_settings

#endif
