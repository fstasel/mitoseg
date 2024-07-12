/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "detection.h"

template <typename T> T quick_select(T arr[], int n) {
    int low, high;
    int median;
    int middle, ll, hh;

    low = 0;
    high = n - 1;
    median = (low + high) / 2;
    while (true) {
        if (high <= low) /* One element only */
            return arr[median];

        if (high == low + 1) { /* Two elements only */
            if (arr[low] > arr[high])
                swap(arr[low], arr[high]);
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position
         * low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            swap(arr[middle], arr[high]);
        if (arr[low] > arr[high])
            swap(arr[low], arr[high]);
        if (arr[middle] > arr[low])
            swap(arr[middle], arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        swap(arr[middle], arr[low + 1]);

        /* Nibble from each end towards middle, swapping items when stuck
         */
        ll = low + 1;
        hh = high;
        while (true) {
            do
                ll++;
            while (arr[low] > arr[ll]);
            do
                hh--;
            while (arr[hh] > arr[low]);

            if (hh < ll)
                break;

            swap(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        swap(arr[low], arr[hh]);

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
    // for no warning
    return arr[0];
}

int isValidMitoSlice_25d(snake25d &s, int slice, dataPacket &dp) {
    int z1, z2;
    int sl, sl_ces_lo, sl_ces_hi, sl_grad_ces_lo, sl_cps_hi;

    // Find z-intersection of data
    z1 = (dp.ces_lo.start_z > dp.ces_hi.start_z) ? dp.ces_lo.start_z
                                                 : dp.ces_hi.start_z;
    z1 = (z1 > dp.grad_ces_lo.start_z) ? z1 : dp.grad_ces_lo.start_z;
    z1 = (z1 > dp.cps_hi.start_z) ? z1 : dp.cps_hi.start_z;
    z1 = (z1 > s.start_z) ? z1 : s.start_z;

    z2 =
        (dp.ces_lo.end_z < dp.ces_hi.end_z) ? dp.ces_lo.end_z : dp.ces_hi.end_z;
    z2 = (z2 < dp.grad_ces_lo.end_z) ? z2 : dp.grad_ces_lo.end_z;
    z2 = (z2 < dp.cps_hi.end_z) ? z2 : dp.cps_hi.end_z;
    z2 = (z2 < s.end_z) ? z2 : s.end_z;

    // Check z-coverage
    if (z1 != s.start_z || z2 != s.end_z) {
        cerr << "Datapacket does not match!" << endl;
        return 0;
    }

    // Averaging range
    int h_range = DT25D_MED_RANGE / 2;
    int start = slice - h_range;
    int end = slice + h_range;
    if (start < 0)
        start = 0;
    if (end > sn25d_t - 1)
        end = sn25d_t - 1;
    int range = end - start + 1;

    int _h_range = h_range, _start = start, _end = end, _range = range;

    // Index shift
    sl_ces_lo = z1 - dp.ces_lo.start_z;
    sl_ces_hi = z1 - dp.ces_hi.start_z;
    sl_grad_ces_lo = z1 - dp.grad_ces_lo.start_z;
    sl_cps_hi = z1 - dp.cps_hi.start_z;
    ///////////////////////////////////////////////////////

    double nx, ny, gx, gy, m;
    int ip, in;
    double dx, dy, d;
    double tx, ty, t;
    int i, j, ii, f, k, w;
    double c;

    // Prepare mitochondrial region
    int n[] = {SN_N};
    cv::Mat region[SN25D_T];
    for (sl = start; sl <= end; sl++) {
        region[sl] = cv::Mat::zeros(cv::Size(dp.ces_lo.ces[sl + sl_ces_lo].w,
                                             dp.ces_lo.ces[sl + sl_ces_lo].h),
                                    CV_8UC1);
    }

    cv::Point p[SN25D_T][SN_N];
    cv::Point *pa;

    for (sl = start; sl <= end; sl++) {
        for (i = 0; i < SN_N; i++) {
            p[sl][i].x = cvRound(s.node[sl][i][0]);
            p[sl][i].y = cvRound(s.node[sl][i][1]);
        }
        pa = p[sl];
        cv::fillPoly(region[sl], (const cv::Point **)&pa, n, 1,
                     cv::Scalar(255));
    }
    h_range = 0;
    start = end = slice;
    range = 1;
    // Compute area
    double area[SN25D_T];
    double avgArea, volume = 0;
    for (sl = start; sl <= end; sl++) {
        area[sl] = 0;
        for (i = 0, j = 1; i < SN_N; i++, j = (j + 1) % SN_N) {
            area[sl] += s.node[sl][i][0] * s.node[sl][j][1] -
                        s.node[sl][i][1] * s.node[sl][j][0];
        }
        area[sl] = 0.5 * fabs(area[sl]);
        volume += area[sl];
    }
    avgArea = volume / range;

    // Compute circumference
    double circum[SN25D_T];
    double avgCircum = 0;
    for (sl = start; sl <= end; sl++) {
        circum[sl] = 0;
        for (i = 0, j = 2; i < SN_N; i += 2, j = (j + 2) % SN_N) {
            dx = s.node[sl][i][0] - s.node[sl][j][0];
            dy = s.node[sl][i][1] - s.node[sl][j][1];
            d = sqrt(dx * dx + dy * dy);
            circum[sl] += d;
        }
        avgCircum += circum[sl];
    }
    avgCircum /= range;

    // Compute curvature (scale variant)
    double curvL = 0;
    double curvG = 0;
    for (j = 0; j < SN_N; j += 2) {
        c = 0;
        for (sl = start; sl <= end; sl++) {
            for (i = 0; i < cvRound(SN_N * DT_CURVL_SEG / 2); i += 2) {
                ip = (i + j + 2) % SN_N;
                in = (i + j + SN_N - 2) % SN_N;
                ii = (i + j) % SN_N;
                gx = s.node[sl][ip][0] - s.node[sl][ii][0];
                gy = s.node[sl][ip][1] - s.node[sl][ii][1];
                d = sqrt(gx * gx + gy * gy);
                gx /= d;
                gy /= d;
                dx = s.node[sl][ii][0] - s.node[sl][in][0];
                dy = s.node[sl][ii][1] - s.node[sl][in][1];
                d = sqrt(dx * dx + dy * dy);
                dx /= d;
                dy /= d;
                tx = gx - dx;
                ty = gy - dy;
                t = sqrt(tx * tx + ty * ty);
                c += t / d;
                // Integral (dt/dd)dd
                if (i == 0)
                    curvG += t;
            }
        }
        c /= range;
        if (curvL < c)
            curvL = c;
    }
    curvG /= range;

    // Signature analysis
    double avgMinMinorRatio;
    int avgMaxNumCriticals = 0;
    double avgMajorAxis = 0, avgMinAxis = 0;
    double avgMinorAxis = 0;
    for (sl = start; sl <= end; sl++) {
        double sign[SN_N][SN_N], sign_smooth[SN_N][SN_N];
        double critical[SN_N][SN_N];
        int nc[SN_N];
        double majorAxis = 0, minAxis = 999999999.;
        int ma1 = 0, ma2 = 0;
        int maxNumCriticals = 0;
        double sign_th = DT_SIGN_RATIO * (circum[sl] / SN_N);
        int window = cvRound(DT_SIGN_SMOOTH * SN_N / 2);
        for (j = 0; j < SN_N; j++) {
            // Compute signature w.r.t. point j
            for (i = 0; i < SN_N; i++) {
                k = (i + j) % SN_N;
                dx = s.node[sl][k][0] - s.node[sl][j][0];
                dy = s.node[sl][k][1] - s.node[sl][j][1];
                sign[j][i] = sqrt(dx * dx + dy * dy);
            }
            // Smooth signature function
            for (i = 0; i < SN_N; i++) {
                c = 0;
                for (w = -window; w <= window; w++) {
                    k = (i + w + SN_N) % SN_N;
                    c += sign[j][k];
                }
                sign_smooth[j][i] = c / (2 * window + 1);
            }
            // Find critical points
            c = critical[j][0] = sign_smooth[j][0];
            nc[j] = 1;
            f = 1;
            for (i = 1; i < SN_N; i++) {
                if ((f && sign_smooth[j][i] < c) ||
                    (!f && sign_smooth[j][i] > c)) {
                    // Store critical pt. if not on border
                    if (!(s.node[sl][j][0] <= 1 || s.node[sl][j][1] <= 1 ||
                          s.node[sl][j][0] >=
                              dp.ces_lo.ces[sl + sl_ces_lo].w - 2 ||
                          s.node[sl][j][1] >=
                              dp.ces_lo.ces[sl + sl_ces_lo].h - 2 ||
                          s.node[sl][i][0] <= 1 || s.node[sl][i][1] <= 1 ||
                          s.node[sl][i][0] >=
                              dp.ces_lo.ces[sl + sl_ces_lo].w - 2 ||
                          s.node[sl][i][1] >=
                              dp.ces_lo.ces[sl + sl_ces_lo].h - 2)) {
                        critical[j][nc[j]] = c;
                        nc[j]++;
                    }
                    // Extract global signature features
                    if (minAxis > c)
                        minAxis = c;
                    if (majorAxis < c) {
                        majorAxis = c;
                        ma1 = j;
                        ma2 = (j + i) % SN_N;
                    }
                    //
                    f = 1 - f;
                }
                c = sign_smooth[j][i];
            }
        }
        for (j = 0; j < SN_N; j++) {
            ii = 0;
            for (i = 0; i < nc[j]; i++) {
                k = (i + 1) % nc[j];
                d = fabs(critical[j][i] - critical[j][k]);
                if (d > sign_th) {
                    ii++;
                }
            }
            if (maxNumCriticals < ii)
                maxNumCriticals = ii;
        }

        // Compute minor axis
        double m1 = 0, m2 = 0;
        double minorAxis;
        dx = s.node[sl][ma1][1] - s.node[sl][ma2][1];
        dy = s.node[sl][ma2][0] - s.node[sl][ma1][0];
        d = sqrt(dx * dx + dy * dy);
        dx /= d;
        dy /= d;
        for (i = 0; i < SN_N; i++) {
            gx = s.node[sl][i][0] - s.node[sl][ma1][0];
            gy = s.node[sl][i][1] - s.node[sl][ma1][1];
            d = gx * dx + gy * dy;
            if (m1 > d)
                m1 = d;
            if (m2 < d)
                m2 = d;
        }
        minorAxis = fabs(m2 - m1);

        // Averaging
        avgMaxNumCriticals += maxNumCriticals;
        avgMajorAxis += majorAxis;
        avgMinAxis += minAxis;
        avgMinorAxis += minorAxis;
    }
    avgMaxNumCriticals = cvRound((double)avgMaxNumCriticals / range);
    avgMajorAxis /= range;
    avgMinAxis /= range;
    avgMinorAxis /= range;
    avgMinMinorRatio = avgMinAxis / avgMinorAxis;
    h_range = _h_range;
    start = _start;
    end = _end;
    range = _range;
    // Compute boundary energy
    double Eboundary[SN25D_T][SN_N];
    double totalEboundary = 0;
    double avgEboundary;
    for (sl = start; sl <= end; sl++) {
        for (i = 0, ip = SN_N - 1, in = 1; i < SN_N;
             i++, ip = (ip + 1) % SN_N, in = (in + 1) % SN_N) {
            nx = s.node[sl][ip][1] - s.node[sl][in][1];
            ny = s.node[sl][in][0] - s.node[sl][ip][0];
            m = sqrt(nx * nx + ny * ny);
            nx /= m;
            ny /= m;
            j = p[sl][i].y * dp.grad_ces_lo.grad_ces[sl + sl_grad_ces_lo].w +
                p[sl][i].x;
            gx = dp.grad_ces_lo.grad_ces[sl + sl_grad_ces_lo].grad_ce[j].x;
            gy = dp.grad_ces_lo.grad_ces[sl + sl_grad_ces_lo].grad_ce[j].y;
            m = sqrt(gx * gx + gy * gy);
            Eboundary[sl][i] = 0;
            if (m > EPS) {
                j = p[sl][i].y * dp.ces_lo.ces[sl + sl_ces_lo].w + p[sl][i].x;
                gx *= dp.ces_lo.ces[sl + sl_ces_lo].ce[j] / m;
                gy *= dp.ces_lo.ces[sl + sl_ces_lo].ce[j] / m;
                Eboundary[sl][i] = fabs(gx * nx + gy * ny);
                totalEboundary += Eboundary[sl][i];
            }
        }
    }
    avgEboundary = totalEboundary / (SN_N * range);
    h_range = 0;
    start = end = slice;
    range = 1;
    // Compute region energy
    double totalEregion = 0;
    double avgEregion;
    for (sl = start; sl <= end; sl++) {
        double xscale = (double)dp.ces_lo.ces[sl + sl_ces_lo].w /
                        dp.ces_hi.ces[sl + sl_ces_hi].w;
        double yscale = (double)dp.ces_lo.ces[sl + sl_ces_lo].h /
                        dp.ces_hi.ces[sl + sl_ces_hi].h;
        int cx, cy;
        for (i = 0; i < (int)dp.cps_hi.cps[sl + sl_cps_hi].size(); i++) {
            cx = cvRound(dp.cps_hi.cps[sl + sl_cps_hi][i].x * xscale);
            cy = cvRound(dp.cps_hi.cps[sl + sl_cps_hi][i].y * yscale);
            if ((region[sl]).at<uchar>(cy, cx)) {
                int j = dp.cps_hi.cps[sl + sl_cps_hi][i].y *
                            dp.ces_hi.ces[sl + sl_ces_hi].w +
                        dp.cps_hi.cps[sl + sl_cps_hi][i].x;
                totalEregion += dp.ces_hi.ces[sl + sl_ces_hi].ce[j];
            }
        }
    }
    avgEregion = totalEregion / volume;
    h_range = _h_range;
    start = _start;
    end = _end;
    range = _range;
    // Compute gap statistics
    double Eavg_boundary[SN_N];
    double eb[SN25D_T];
    for (i = 0; i < SN_N; i++) {
        for (sl = start; sl <= end; sl++) {
            eb[sl] = Eboundary[sl][i];
        }
        Eavg_boundary[i] = quick_select(&eb[start], range);
    }
    double maxGap = 0;
    double totalGap = 0;
    double gapRatio;
    double borderGapRatio = 0;
    int countGap = 0, countBorderGap = 0;
    bool gapOnBorder = false;
    int borderGaps[SN25D_T][4];
    int sumBG;
    for (sl = start; sl <= end; sl++)
        borderGaps[sl][0] = borderGaps[sl][1] = borderGaps[sl][2] =
            borderGaps[sl][3] = 0;
    f = 0;
    c = 0;
    i = 0;
    j = SN_N + 1;
    while (i < j) {
        ii = i % SN_N;
        if (Eavg_boundary[ii] < DT_GAP_T) {
            if (!f) {
                j++;
            } else {
                // Count gaps on borders
                for (sl = start; sl <= end; sl++) {
                    if (s.node[sl][ii][0] <= 1)
                        borderGaps[sl][0] = 1;
                    if (s.node[sl][ii][1] <= 1)
                        borderGaps[sl][1] = 1;
                    if (s.node[sl][ii][0] >=
                        dp.ces_lo.ces[sl + sl_ces_lo].w - 2)
                        borderGaps[sl][2] = 1;
                    if (s.node[sl][ii][1] >=
                        dp.ces_lo.ces[sl + sl_ces_lo].h - 2)
                        borderGaps[sl][3] = 1;
                }
                // gap on border?
                int gapOnBorderCount = 0;
                for (sl = start; sl <= end; sl++) {
                    if (s.node[sl][ii][0] <= 1 || s.node[sl][ii][1] <= 1 ||
                        s.node[sl][ii][0] >=
                            dp.ces_lo.ces[sl + sl_ces_lo].w - 2 ||
                        s.node[sl][ii][1] >=
                            dp.ces_lo.ces[sl + sl_ces_lo].h - 2)
                        gapOnBorderCount++;
                }
                if (gapOnBorderCount > range / 2)
                    gapOnBorder = true;
                //
                d = 0;
                for (sl = start; sl <= end; sl++) {
                    ip = (i + SN_N - 1) % SN_N;
                    in = (i + 1) % SN_N;
                    dx = s.node[sl][ii][0] - s.node[sl][ip][0];
                    dy = s.node[sl][ii][1] - s.node[sl][ip][1];
                    d += sqrt(dx * dx + dy * dy);
                    dx = s.node[sl][ii][0] - s.node[sl][in][0];
                    dy = s.node[sl][ii][1] - s.node[sl][in][1];
                    d += sqrt(dx * dx + dy * dy);
                }
                d *= 0.5;
                d /= range;
                totalGap += d;
                c += d;
            }
        } else {
            if (gapOnBorder) {
                borderGapRatio += c;
                gapOnBorder = false;
            }
            if (maxGap < c)
                maxGap = c;
            if (c > 0) {
                countGap++;
                c = 0;
            }
            f = 1;
        }
        i++;
        if (i == SN_N && !f) {
            maxGap = totalGap = avgCircum;
            countGap = 1;
            sumBG = 0;
            for (sl = start; sl <= end; sl++)
                sumBG += borderGaps[sl][0] + borderGaps[sl][1] +
                         borderGaps[sl][2] + borderGaps[sl][3];
            borderGapRatio = (sumBG > 0) ? avgCircum : 0;
            break;
        }
    }
    borderGapRatio /= avgCircum;
    gapRatio = totalGap / avgCircum;
    countBorderGap = 0;
    for (j = 0; j < 4; j++) {
        sumBG = 0;
        for (sl = start; sl <= end; sl++)
            sumBG += borderGaps[sl][j];
        if (sumBG > range / 2)
            countBorderGap++;
    }

    // Validator
    bool c1 = (avgEboundary >= DT_BOUNDARY_T);
    bool c2 = (avgEregion >= DT_REGION_T);
    bool c3 = (maxGap <= DT_GAPMAX_T);
    bool c4 = (totalGap <= DT_GAPTOTAL_T);
    bool c5 = (gapRatio <= DT_GAPTOTALRATIO_T);
    bool c6 = (borderGapRatio <= DT_GAPBORDERRATIO_T);
    bool c7 = (countBorderGap <= DT_GAPBORDERCOUNT_T);
    bool c8 = (curvL <= DT_CURVL_T);
    bool c9 = (curvG <= DT_CURVG_T);
    bool c10 = (avgMaxNumCriticals <= DT_SIGN_MAXNUMCRIT_T);
    bool c11 = (avgMajorAxis <= DT_MAJORAXIS_T);
    bool c16 = (avgMinorAxis >= DT_MINORAXIS_T);
    bool c12 = (avgMinAxis >= DT_MINAXIS_T);
    bool c13 = (avgMinMinorRatio >= DT_MINMINOR_RATIO_T);
    bool c14 = (avgArea >= DT_MINAREA);
    bool c15 = (avgArea <= DT_MAXAREA);
    int is_Valid = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9 && c10 &&
                   c11 && c12 && c13 && c14 && c15 && c16;

    // Report
    if (DT_REPORT) {
        cout << "------------" << endl;
        cout << boost::format("Eboundary = %f %c") % avgEboundary %
                    (c1 ? ' ' : '-')
             << endl;
        cout << boost::format("Eregion = %f %c") % avgEregion % (c2 ? ' ' : '-')
             << endl;
        cout << boost::format("Maxgap = %f (%f nm) %c") % maxGap %
                    (maxGap * TFACTOR * LE_SSIZE_LO) % (c3 ? ' ' : '-')
             << endl;
        cout << boost::format("Totalgap = %f (%f nm) %c") % totalGap %
                    (totalGap * TFACTOR * LE_SSIZE_LO) % (c4 ? ' ' : '-')
             << endl;
        cout << boost::format("Gap ratio = %f %c") % gapRatio % (c5 ? ' ' : '-')
             << endl;
        cout << boost::format("Bordergap ratio = %f %c") % borderGapRatio %
                    (c6 ? ' ' : '-')
             << endl;
        cout << boost::format("Gap count = %d") % countGap << endl;
        cout << boost::format("Bordergap count = %d %c") % countBorderGap %
                    (c7 ? ' ' : '-')
             << endl;
        cout << boost::format("Curvature Local = %f %c") % curvL %
                    (c8 ? ' ' : '-')
             << endl;
        cout << boost::format("Curvature Global = %f %c") % curvG %
                    (c9 ? ' ' : '-')
             << endl;
        cout << boost::format("Area = %f (%f nm2) %c") % avgArea %
                    (avgArea * (TFACTOR * TFACTOR) *
                     (LE_SSIZE_LO * LE_SSIZE_LO)) %
                    ((c14 || c15) ? ' ' : '-')
             << endl;
        cout << boost::format("Circumference = %f (%f nm)") % avgCircum %
                    (avgCircum * TFACTOR * LE_SSIZE_LO)
             << endl;
        cout << boost::format("Major axis = %f (%f nm) %c") % avgMajorAxis %
                    (avgMajorAxis * TFACTOR * LE_SSIZE_LO) % (c11 ? ' ' : '-')
             << endl;
        cout << boost::format("Minor axis = %f (%f nm) %c") % avgMinorAxis %
                    (avgMinorAxis * TFACTOR * LE_SSIZE_LO) % (c16 ? ' ' : '-')
             << endl;
        cout << boost::format("Min. axis = %f (%f nm) %c") % avgMinAxis %
                    (avgMinAxis * TFACTOR * LE_SSIZE_LO) % (c12 ? ' ' : '-')
             << endl;
        cout << boost::format("Min/Minor ratio = %f %c") % avgMinMinorRatio %
                    (c13 ? ' ' : '-')
             << endl;
        cout << boost::format("Max nc = %d %c") % avgMaxNumCriticals %
                    (c10 ? ' ' : '-')
             << endl;
        cout << "------------" << endl;
    }

    return is_Valid;
}
