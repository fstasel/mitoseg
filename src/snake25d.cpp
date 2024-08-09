/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "snake25d.h"

#include "detection.h"

sem_t sem2;
pthread_mutex_t addlock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t tnslock = PTHREAD_MUTEX_INITIALIZER;

curveStack loadCurveStack(int start_z, int end_z, const string &tag) {
    int s, i;
    curveStack cstack;
    cstack.start_z = start_z;
    cstack.end_z = end_z;
    cstack.t = end_z - start_z + 1;

    for (i = 0, s = start_z; s <= end_z; i++, s++) {
        cstack.cstack.push_back(loadCurveStructure(s, tag));
    }
    return cstack;
}

curvePoints<double> createCurvePoints(curveList &c) {
    uint i;
    curvePoints<double> cp;

    for (i = 0; i < c.size(); i++) {
        curvePoints<double> t = getCurvePoints(c[i], 1.0);
        cp.insert(cp.end(), t.begin(), t.end());
    }
    return cp;
}

curvePoints<int> createCurvePoints(curveList &c, int w, int h) {
    curvePoints<double> temp = createCurvePoints(c);
    curvePoints<int> cp(temp.size());

    uint j;
    for (j = 0; j < cp.size(); j++) {
        cp[j].x = clamp(cvRound(temp[j].x), 0, w - 1);
        cp[j].y = clamp(cvRound(temp[j].y), 0, h - 1);
    }

    uint msize = w * h;
    bool *map = new bool[msize]();

    for (j = 0; j < cp.size(); j++)
        map[cp[j].y * w + cp[j].x] = true;

    cp.clear();
    for (j = 0; j < msize; j++) {
        if (map[j]) {
            cp.push_back({(int)j % w, (int)j / w});
        }
    }
    delete[] map;

    return cp;
}

curvePointsStack<double> createCurvePointsStack(curveStack &cstack) {
    curvePointsStack<double> cps;
    cps.cps.resize(cstack.t);
    cps.start_z = cstack.start_z;
    cps.end_z = cstack.end_z;
    cps.t = cstack.t;
    for (int i = 0; i < cps.t; i++) {
        cps.cps[i] = createCurvePoints(cstack.cstack[i]);
    }
    return cps;
}

curvePointsStack<int> createCurvePointsStack(curveStack &cstack, int w, int h) {
    curvePointsStack<int> cps;
    cps.cps.resize(cstack.t);
    cps.start_z = cstack.start_z;
    cps.end_z = cstack.end_z;
    cps.t = cstack.t;
    for (int i = 0; i < cps.t; i++) {
        cps.cps[i] = createCurvePoints(cstack.cstack[i], w, h);
    }
    return cps;
}

curveEnergyMap createCurveEnergyMap(curveList &cl, const int w, const int h) {
    curveEnergyMap ce;
    ce.w = w;
    ce.h = h;
    ce.ce = getCurveEnergyImage(cl, w, h);
    return ce;
}

curveEnergyStack createCurveEnergyStack(curveStack &cs, const int w,
                                        const int h) {
    curveEnergyStack ces;
    ces.ces.resize(cs.t);
    ces.start_z = cs.start_z;
    ces.end_z = cs.end_z;
    ces.t = cs.t;

    for (int i = 0; i < ces.t; i++) {
        ces.ces[i] = createCurveEnergyMap(cs.cstack[i], w, h);
    }
    return ces;
}

curveEnergyGradientMap createCurveEnergyGradientMap(curveEnergyMap &ce) {
    curveEnergyGradientMap grad_ce;
    grad_ce.grad_ce = getEnergyGradient(ce.ce, ce.w, ce.h);
    grad_ce.w = ce.w;
    grad_ce.h = ce.h;
    return grad_ce;
}

curveEnergyGradientStack createCurveEnergyGradientStack(curveEnergyStack &ces) {
    curveEnergyGradientStack grad_ces;
    grad_ces.grad_ces.resize(ces.t);
    grad_ces.start_z = ces.start_z;
    grad_ces.end_z = ces.end_z;
    grad_ces.t = ces.t;

    for (int i = 0; i < ces.t; i++) {
        grad_ces.grad_ces[i] = createCurveEnergyGradientMap(ces.ces[i]);
    }
    return grad_ces;
}

void getMapSize(cv::Size ridgeImageSize, const int ssize, int *w, int *h) {
    *w = (ridgeImageSize.width / ssize) +
         ((ridgeImageSize.width % ssize) ? 1 : 0);
    *h = (ridgeImageSize.height / ssize) +
         ((ridgeImageSize.height % ssize) ? 1 : 0);
}

void compSnakeSliceArea(snake25d &s, const int z) {
    int k, j;
    double area = 0;
    for (k = 0, j = 1; k < SN_N; k++, j = (j + 1) % SN_N) {
        area += s.node[z][k][0] * s.node[z][j][1] -
                s.node[z][k][1] * s.node[z][j][0];
    }
    s.area[z] = 0.5 * area;
    if (s.area[z] < 0)
        s.area[z] = -s.area[z];
}

void compSnakeArea(snake25d &s) {
    for (int z = 0; z < s.t; z++) {
        compSnakeSliceArea(s, z);
    }
}

void compApproxIntersectionArea(double p1[][2], int n1, double p2[][2], int n2,
                                double *a1, double *a2, double *is,
                                double precision) {
    int i, x1, y1, x2, y2, w, h;
    cv::Point *k1 = nullptr;
    cv::Point *k2 = nullptr;
    if (n1)
        k1 = new cv::Point[n1];
    if (n2)
        k2 = new cv::Point[n2];
    for (i = 0; i < n1; i++) {
        k1[i].x = (int)(p1[i][0] * precision);
        k1[i].y = (int)(p1[i][1] * precision);
    }
    for (i = 0; i < n2; i++) {
        k2[i].x = (int)(p2[i][0] * precision);
        k2[i].y = (int)(p2[i][1] * precision);
    }
    if (n1) {
        x1 = x2 = k1[0].x;
        y1 = y2 = k1[0].y;
    } else {
        x1 = x2 = k2[0].x;
        y1 = y2 = k2[0].y;
    }
    for (i = 0; i < n1; i++) {
        x1 = (x1 > k1[i].x) ? k1[i].x : x1;
        y1 = (y1 > k1[i].y) ? k1[i].y : y1;
        x2 = (x2 < k1[i].x) ? k1[i].x : x2;
        y2 = (y2 < k1[i].y) ? k1[i].y : y2;
    }
    for (i = 0; i < n2; i++) {
        x1 = (x1 > k2[i].x) ? k2[i].x : x1;
        y1 = (y1 > k2[i].y) ? k2[i].y : y1;
        x2 = (x2 < k2[i].x) ? k2[i].x : x2;
        y2 = (y2 < k2[i].y) ? k2[i].y : y2;
    }
    w = x2 - x1 + 3;
    h = y2 - y1 + 3;
    for (i = 0; i < n1; i++) {
        k1[i].x -= (x1 - 1);
        k1[i].y -= (y1 - 1);
    }
    for (i = 0; i < n2; i++) {
        k2[i].x -= (x1 - 1);
        k2[i].y -= (y1 - 1);
    }
    cv::Mat buf1 = cv::Mat::zeros(h, w, CV_8UC1);
    cv::Mat buf2 = cv::Mat::zeros(h, w, CV_8UC1);
    if (n1)
        cv::fillPoly(buf1, (const cv::Point **)&k1, &n1, 1,
                     cv::Scalar(255, 255, 255));
    if (n2)
        cv::fillPoly(buf2, (const cv::Point **)&k2, &n2, 1,
                     cv::Scalar(255, 255, 255));
    int c1 = cv::countNonZero(buf1);
    int c2 = cv::countNonZero(buf2);
    cv::bitwise_and(buf1, buf2, buf1);
    int inter = cv::countNonZero(buf1);
    *a1 = c1 / (precision * precision);
    *a2 = c2 / (precision * precision);
    *is = inter / (precision * precision);
    if (k1)
        delete[] k1;
    if (k2)
        delete[] k2;
}

double getDiceSimilarity2d(double p1[][2], int n1, double p2[][2], int n2,
                           double precision) {
    double a1, a2, is, ds;
    compApproxIntersectionArea(p1, n1, p2, n2, &a1, &a2, &is, precision);
    ds = 2.0 * is / (a1 + a2);
    return ds;
}

double getDiceSimilaritySnake25d(snake25d &s1, snake25d &s2, double precision) {
    int z1 = (s1.start_z < s2.start_z) ? s1.start_z : s2.start_z;
    int z2 = (s1.end_z > s2.end_z) ? s1.end_z : s2.end_z;
    int z;
    double a1, a2, is;
    double v1 = 0, v2 = 0, iv = 0, ds = 1.0;

    for (z = z1; z <= z2; z++) {
        if (z >= s1.start_z && z >= s2.start_z && z <= s1.end_z &&
            z <= s2.end_z) {
            compApproxIntersectionArea(s1.node[z - s1.start_z], SN_N,
                                       s2.node[z - s2.start_z], SN_N, &a1, &a2,
                                       &is, precision);
        } else if (z >= s1.start_z && z <= s1.end_z) {
            compApproxIntersectionArea(s1.node[z - s1.start_z], SN_N, nullptr,
                                       0, &a1, &a2, &is, precision);
        } else {
            compApproxIntersectionArea(nullptr, 0, s2.node[z - s2.start_z],
                                       SN_N, &a1, &a2, &is, precision);
        }
        v1 += a1;
        v2 += a2;
        iv += is;
    }
    ds = 2.0 * iv / (v1 + v2);
    return ds;
}

snake25d initSnake25d(const double x, const double y, const int start_z,
                      const int end_z, const double init_r, const int w,
                      const int h, const double ar_z) {
    int i, z;
    double rad, sx, sy;

    snake25d s;
    s.start_z = start_z;
    s.end_z = end_z;
    s.t = end_z - start_z + 1;
    s.aspectRatio_z = ar_z;
    s.status = SNAKE_GROWING;
    s.validity = 0;

    for (i = 0; i < SN_N; i++) {
        rad = 2 * PI * i / SN_N;
        sx = x + init_r * cos(rad);
        sy = y + init_r * sin(rad);
        // Viewpoint constraints & correction
        if (sx < 0)
            sx = 0;
        if (sy < 0)
            sy = 0;
        if (sx >= w)
            sx = w - 1;
        if (sy >= h)
            sy = h - 1;

        for (z = 0; z < s.t; z++) {
            s.node[z][i][0] = sx;
            s.node[z][i][1] = sy;
            //
            s.grad_Eint[z][i][0] = 0;
            s.grad_Eint[z][i][1] = 0;
            s.grad_Eext[z][i][0] = 0;
            s.grad_Eext[z][i][1] = 0;
            s.d1_xy[z][i][0] = 0;
            s.d1_xy[z][i][1] = 0;
            s.d2_xy[z][i][0] = 0;
            s.d2_xy[z][i][1] = 0;
            s.d4_xy[z][i][0] = 0;
            s.d4_xy[z][i][1] = 0;
            s.d1_z[z][i][0] = 0;
            s.d1_z[z][i][1] = 0;
            s.d2_z[z][i][0] = 0;
            s.d2_z[z][i][1] = 0;
            s.d4_z[z][i][0] = 0;
            s.d4_z[z][i][1] = 0;
            //
            s.isValid[z] = SNAKE_INVALID;
        }
    }

    // Init area
    compSnakeArea(s);
    return s;
}

void compSnakeDerivatives(snake25d &s) {
    int i, t;
    int k1, k2, k3, k4, k5;
    int t1, t2, t3, t4, t5;

    for (t = 0; t < s.t; t++) {
        // z indices
        t1 = t - 2;
        t2 = t - 1;
        t3 = t;
        t4 = t + 1;
        t5 = t + 2;
        if (t1 < 0)
            t1 = 0;
        if (t2 < 0)
            t2 = 0;
        if (t4 >= s.t)
            t4 = s.t - 1;
        if (t5 >= s.t)
            t5 = s.t - 1;

        // init xy indices
        k1 = SN_N - 2;
        k2 = SN_N - 1;
        k3 = 0;
        k4 = 1;
        k5 = 2;
        for (i = 0; i < SN_N; i++) {
            // xy derivatives
            s.d1_xy[t][i][0] = 0.5 * (s.node[t][k2][0] - s.node[t][k4][0]);
            s.d1_xy[t][i][1] = 0.5 * (s.node[t][k2][1] - s.node[t][k4][1]);
            s.d2_xy[t][i][0] = 0.25 * (-s.node[t][k2][0] +
                                       2 * s.node[t][k3][0] - s.node[t][k4][0]);
            s.d2_xy[t][i][1] = 0.25 * (-s.node[t][k2][1] +
                                       2 * s.node[t][k3][1] - s.node[t][k4][1]);
            s.d4_xy[t][i][0] =
                0.0625 * (s.node[t][k1][0] - 4 * s.node[t][k2][0] +
                          6 * s.node[t][k3][0] - 4 * s.node[t][k4][0] +
                          s.node[t][k5][0]);
            s.d4_xy[t][i][1] =
                0.0625 * (s.node[t][k1][1] - 4 * s.node[t][k2][1] +
                          6 * s.node[t][k3][1] - 4 * s.node[t][k4][1] +
                          s.node[t][k5][1]);

            // z derivatives
            s.d1_z[t][i][0] =
                s.aspectRatio_z * 0.5 * (s.node[t2][i][0] - s.node[t4][i][0]);
            s.d1_z[t][i][1] =
                s.aspectRatio_z * 0.5 * (s.node[t2][i][1] - s.node[t4][i][1]);
            s.d2_z[t][i][0] =
                s.aspectRatio_z * 0.25 *
                (-s.node[t2][i][0] + 2 * s.node[t3][i][0] - s.node[t4][i][0]);
            s.d2_z[t][i][1] =
                s.aspectRatio_z * 0.25 *
                (-s.node[t2][i][1] + 2 * s.node[t3][i][1] - s.node[t4][i][1]);
            s.d4_z[t][i][0] = s.aspectRatio_z * 0.0625 *
                              (s.node[t1][i][0] - 4 * s.node[t2][i][0] +
                               6 * s.node[t3][i][0] - 4 * s.node[t4][i][0] +
                               s.node[t5][i][0]);
            s.d4_z[t][i][1] = s.aspectRatio_z * 0.0625 *
                              (s.node[t1][i][1] - 4 * s.node[t2][i][1] +
                               6 * s.node[t3][i][1] - 4 * s.node[t4][i][1] +
                               s.node[t5][i][1]);

            // update xy indices
            k1 = (k1 + 1) % SN_N;
            k2 = (k2 + 1) % SN_N;
            k3 = (k3 + 1) % SN_N;
            k4 = (k4 + 1) % SN_N;
            k5 = (k5 + 1) % SN_N;
        }
    }
}

void compIntEnergyGrad(snake25d &s, const double w_tens_xy,
                       const double w_curv_xy, const double w_tens_z,
                       const double w_curv_z) {
    int i, t;
    for (t = 0; t < s.t; t++) {
        for (i = 0; i < SN_N; i++) {
            // Gradient vector of internal energy
            s.grad_Eint[t][i][0] =
                w_tens_xy * s.d2_xy[t][i][0] + w_curv_xy * s.d4_xy[t][i][0] +
                w_tens_z * s.d2_z[t][i][0] + w_curv_z * s.d4_z[t][i][0];

            s.grad_Eint[t][i][1] =
                w_tens_xy * s.d2_xy[t][i][1] + w_curv_xy * s.d4_xy[t][i][1] +
                w_tens_z * s.d2_z[t][i][1] + w_curv_z * s.d4_z[t][i][1];
        }
    }
}

void compExtEnergyGrad(snake25d &s, curveEnergyGradientStack &grad_ces,
                       const double w_curve, const double w_inf) {
    int i, j, x, y, t, k;
    double m, nx, ny;

    for (t = 0, k = s.start_z - grad_ces.start_z; t < s.t; t++, k++) {
        for (i = 0; i < SN_N; i++) {
            x = cvRound(s.node[t][i][0]);
            y = cvRound(s.node[t][i][1]);
            j = y * grad_ces.grad_ces[k].w + x;
            m = sqrt(s.d1_xy[t][i][0] * s.d1_xy[t][i][0] +
                     s.d1_xy[t][i][1] * s.d1_xy[t][i][1]);
            nx = -s.d1_xy[t][i][1] / m;
            ny = s.d1_xy[t][i][0] / m;
            // Compute the gradient of external energy
            s.grad_Eext[t][i][0] =
                -w_curve * grad_ces.grad_ces[k].grad_ce[j].x - w_inf * nx;
            s.grad_Eext[t][i][1] =
                -w_curve * grad_ces.grad_ces[k].grad_ce[j].y - w_inf * ny;
        }
    }
}

void equidistantCorrectionUV(snake25d &s, double uv[][SN_N][2],
                             double uv_corrected[][SN_N][2]) {
    int t, j, jp, jn;
    double dx, dy, dm, cx, cy, c;

    for (t = 0; t < s.t; t++) {
        for (j = 0, jp = SN_N - 1, jn = 1; j < SN_N;
             j++, jp = (jp + 1) % SN_N, jn = (jn + 1) % SN_N) {
            dx = s.node[t][jn][0] - s.node[t][jp][0];
            dy = s.node[t][jn][1] - s.node[t][jp][1];
            dm = dx * dx + dy * dy;
            if (dm < EPS) {
                c = 0;
            } else {
                cx = 0.5 * (s.node[t][jn][0] + s.node[t][jp][0]);
                cy = 0.5 * (s.node[t][jn][1] + s.node[t][jp][1]);
                c = (dx * (cx - (s.node[t][j][0] + uv[t][j][0])) +
                     dy * (cy - (s.node[t][j][1] + uv[t][j][1]))) /
                    dm;
            }
            uv_corrected[t][j][0] = uv[t][j][0] + c * dx;
            uv_corrected[t][j][1] = uv[t][j][1] + c * dy;
        }
    }
}

void compUpdateVectors(snake25d &s, double uv[][SN_N][2]) {
    int t, j;
    double ux, uy, m, mx, gamma;
    for (t = 0; t < s.t; t++) {
        mx = 0;
        for (j = 0; j < SN_N; j++) {
            uv[t][j][0] = ux = s.grad_Eint[t][j][0] + s.grad_Eext[t][j][0];
            uv[t][j][1] = uy = s.grad_Eint[t][j][1] + s.grad_Eext[t][j][1];
            m = ux * ux + uy * uy;
            if (mx < m)
                mx = m;
        }
        mx = sqrt(mx);
        // Compute time discretization parameter gamma
        if (mx > EPS)
            gamma = 1.0 / mx;
        else
            gamma = 0;
        for (j = 0; j < SN_N; j++) {
            uv[t][j][0] *= -SN25D_K * gamma;
            uv[t][j][1] *= -SN25D_K * gamma;
            m = cv::sqrt(uv[t][j][0] * uv[t][j][0] + uv[t][j][1] * uv[t][j][1]);
            m = clamp(m, 0.0, 1.0);
            double ac = alphaCoefficient(m);
            uv[t][j][0] *= (ac / m);
            uv[t][j][1] *= (ac / m);
        }
    }
}

double alphaCoefficient(double x) {
    return cv::pow(1.0 - cv::pow(1.0 - x, 1.0 / (1.0 - SN25D_ALPHA)),
                   1.0 - SN25D_ALPHA);
}

void updateSnake(snake25d &s, double uv[][SN_N][2], const int w_1,
                 const int h_1, double mov[][SN_N][2], double *maxMov) {
    int t, i;
    double px, py, dx, dy, m;
    for (t = 0; t < s.t; t++) {
        for (i = 0; i < SN_N; i++) {
            px = s.node[t][i][0];
            py = s.node[t][i][1];
            // Update snake
            s.node[t][i][0] += uv[t][i][0];
            s.node[t][i][1] += uv[t][i][1];
            // Viewpoint constraints & correction
            s.node[t][i][0] = clamp(s.node[t][i][0], 0.0, (double)w_1);
            s.node[t][i][1] = clamp(s.node[t][i][1], 0.0, (double)h_1);
            // Update movement
            dx = (mov[t][i][0] += s.node[t][i][0] - px);
            dy = (mov[t][i][1] += s.node[t][i][1] - py);
            m = dx * dx + dy * dy;
            if (*maxMov < m)
                *maxMov = m;
        }
    }
}

void resetMovement(double movement[][SN_N][2], double *maxMovement) {
    fill(&movement[0][0][0], &movement[0][0][0] + SN25D_T * SN_N * 2, 0);
    *maxMovement = 0;
}

snake25d startSnake25d(curveEnergyGradientStack &grad_ces,
                       const double w_tens_xy, const double w_curv_xy,
                       const double w_tens_z, const double w_curv_z,
                       const double w_curve, const double w_inf,
                       const double aspectRatio_z, snake25d *initialSnake,
                       const double initx, const double inity,
                       const int start_z, int end_z) {
    snake25d snake;
    double uv[SN25D_T][SN_N][2];
    double uv_corrected[SN25D_T][SN_N][2];
    double movement[SN25D_T][SN_N][2], maxMovement;

    int i, z;

    int w = grad_ces.grad_ces[0].w - 1;
    int h = grad_ces.grad_ces[0].h - 1;
    int w_1 = w - 1;
    int h_1 = h - 1;

    // Initiate snake
    snake = initialSnake ? *initialSnake
                         : initSnake25d(initx, inity, start_z, end_z,
                                        SN25D_INITR, w, h, aspectRatio_z);
    snake.status = SNAKE_GROWING;
    // Initialize iteration counter
    i = 0;
    // Reset cumulative snake movement
    resetMovement(movement, &maxMovement);
    // Main loop
    while (true) {
        // Compute Snake Derivatives
        compSnakeDerivatives(snake);
        // Compute Internal Energy
        compIntEnergyGrad(snake, w_tens_xy, w_curv_xy, w_tens_z, w_curv_z);
        // Compute External Energy
        compExtEnergyGrad(snake, grad_ces, w_curve, w_inf);
        // Compute Update Vectors
        compUpdateVectors(snake, uv);
        // Equidistant correction
        equidistantCorrectionUV(snake, uv, uv_corrected);
        // Update snake
        updateSnake(snake, uv_corrected, w_1, h_1, movement, &maxMovement);
        // Compute area
        compSnakeArea(snake);
        // -----------------
        // Stopping criteria
        // -----------------
        // Max iter.
        if (i > SN25D_MAXITER) {
            snake.status = SNAKE_MAXITER;
            break;
        }
        // Min. area constraint
        for (z = 0; z < snake.t; z++) {
            if (snake.area[z] < SN25D_MINAREA) {
                snake.status = SNAKE_MINAREA;
                break;
            }
        }
        if (snake.status == SNAKE_MINAREA) {
            break;
        }
        // Max. area constraint
        for (z = 0; z < snake.t; z++) {
            if (snake.area[z] > SN25D_MAXAREA) {
                snake.status = SNAKE_MAXAREA;
                break;
            }
        }
        if (snake.status == SNAKE_MAXAREA) {
            break;
        }

        // Convergence criteria
        if (i % SN25D_SHORTTERM_CONV_ITER == SN25D_SHORTTERM_CONV_ITER - 1 ||
            i % SN25D_LONGTERM_CONV_ITER == SN25D_LONGTERM_CONV_ITER - 1) {
            if (maxMovement < SN25D_SHORTTERM_CONV ||
                maxMovement < SN25D_LONGTERM_CONV) {
                snake.status = SNAKE_CONVERGED;
                break;
            }
            resetMovement(movement, &maxMovement);
        }
        // -----------------
        // Iterate
        i++;
    }

    //  Return the output snake
    return snake;
}

int inflationSnake25d(dataPacket &dp, const double initx, const double inity,
                      const int start_z, int end_z, const double w_tens_xy,
                      const double w_curv_xy, const double w_tens_z,
                      const double w_curv_z, const double w_curve,
                      const double w_inf_min, const double w_inf_max,
                      const double w_inf_step, const double aspectRatio_z,
                      snake25dList &outputSnakeList) {
    snake25d snake, outputSnake;
    double w_inf, sim;
    bool flag = false;
    int i = 0;

    for (w_inf = w_inf_min; w_inf <= w_inf_max; w_inf += w_inf_step) {
        outputSnake = startSnake25d(
            dp.grad_ces_lo, w_tens_xy, w_curv_xy, w_tens_z, w_curv_z, w_curve,
            w_inf, aspectRatio_z, nullptr, initx, inity, start_z, end_z);
        if (outputSnake.status == SNAKE_CONVERGED) {
            // Check similarity
            sim = flag ? getDiceSimilaritySnake25d(snake, outputSnake) : 0.0;
            // If we have a different snake
            if (sim <= SN25D_INF_CONV) {
                // Validation
                validateSnake25d(dp, outputSnake);
                // Add to list
                pthread_mutex_lock(&addlock);
                outputSnakeList.push_back(outputSnake);
                pthread_mutex_unlock(&addlock);
                // Report validity
                if (DT_REPORT) {
                    reportValidation(dp, outputSnake);
                }
                i++;
                snake = outputSnake;
                flag = true;
                // Converged --save snake and restart with greater
                // w_inf
            }
        }
        // Over inflation --stop
        if (outputSnake.status == SNAKE_MAXAREA)
            break;
        // Insufficient Inflation --restart with greater w_inf
        if (outputSnake.status == SNAKE_MAXITER) {
            flag = false;
        }
        // Vanished --restart with greater w_inf
        if (outputSnake.status == SNAKE_MINAREA) {
            flag = false;
        }
    }

    return i;
}

void reportValidation(dataPacket &dp, snake25d &s) {
    for (int t = 0; t < s.t; t++) {
        isValidMitoSlice_25d(s, t, dp);
        cv::waitKey();
    }
}

void validateSnake25d(dataPacket &dp, snake25d &s) {
    int v, n_v = 0, n_iv = 0;
    for (int t = 0; t < s.t; t++) {
        v = isValidMitoSlice_25d(s, t, dp);
        if (v) {
            s.isValid[t] = SNAKE_VALID;
            n_v++;
        } else {
            s.isValid[t] = SNAKE_INVALID;
            n_iv++;
        }
    }
    s.validity = (double)n_v / (n_v + n_iv);
}

vector<cv::Point> findInitialPointsWithinZRange(dataPacket &dp, int start_z,
                                                int end_z) {
    int z, zci, zpi, ci, c;
    int i, j, x, y;
    double cx1, cx2, cy1, cy2, midx, midy, dx, dy, dm;
    double px, py;
    int n_avg, n_pts = 0, n_out;

    // Initialize memory
    int msize = 0;
    for (z = start_z; z <= end_z; z++)
        if (z >= dp.cs_lo.start_z && z <= dp.cs_lo.end_z)
            msize += dp.cs_lo.cstack[z - dp.cs_lo.start_z].size();
    double(*rawPoints)[2] = new double[msize][2];
    double(*avgPoints)[2] = new double[msize][2];
    curve **assocCurve = new curve *[msize];
    int *neighbors = new int[msize];
    vector<bool> visit(msize);
    double eps = SN25D_INITPTS_EPS;
    int minPts = SN25D_INITPTS_MIN;
    double distsq = eps * eps;
    double initloc = 2.0 * SN25D_INITR;

    // Determine the possible initial points
    for (z = start_z; z <= end_z; z++) {
        if (z >= dp.cs_lo.start_z && z <= dp.cs_lo.end_z &&
            z >= dp.cps_lo.start_z && z <= dp.cps_lo.end_z) {
            zci = z - dp.cs_lo.start_z;
            zpi = z - dp.cps_lo.start_z;
            ci = 0;
            for (c = 0; c < (int)dp.cs_lo.cstack[zci].size(); c++) {
                // Store initial point
                if (dp.cs_lo.cstack[zci][c].h != 0) {
                    cx1 = dp.cps_lo.cps[zpi][ci].x;
                    cy1 = dp.cps_lo.cps[zpi][ci].y;
                    cx2 =
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len - 1]
                            .x;
                    cy2 =
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len - 1]
                            .y;
                    midx = (cx1 + cx2) / 2.0;
                    midy = (cy1 + cy2) / 2.0;
                    dx =
                        midx -
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len / 2]
                            .x;
                    dy =
                        midy -
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len / 2]
                            .y;
                    dm = sqrt(dx * dx + dy * dy);
                    dx *= initloc / dm;
                    dy *= initloc / dm;
                    rawPoints[n_pts][0] =
                        dx +
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len / 2]
                            .x;
                    rawPoints[n_pts][1] =
                        dy +
                        dp.cps_lo.cps[zpi][ci + dp.cs_lo.cstack[zci][c].len / 2]
                            .y;
                    assocCurve[n_pts] = &dp.cs_lo.cstack[zci][c];
                    n_pts++;
                }

                ci += dp.cs_lo.cstack[zci][c].len;
            }
        }
    }

    // Density-based clustering
    for (i = 0; i < n_pts; i++) {
        c = 0;
        for (j = 0; j < n_pts; j++) {
            dx = rawPoints[i][0] - rawPoints[j][0];
            dy = rawPoints[i][1] - rawPoints[j][1];
            dm = dx * dx + dy * dy;
            if (dm <= distsq &&
                isInside(*assocCurve[i], rawPoints[j][0], rawPoints[j][1]) &&
                isInside(*assocCurve[j], rawPoints[i][0], rawPoints[i][1])) {
                c++;
            }
        }
        neighbors[i] = c;
    }
    for (i = 0; i < n_pts - 1; i++) {
        for (j = i + 1; j < n_pts; j++) {
            if (neighbors[i] < neighbors[j]) {
                swap(neighbors[i], neighbors[j]);
                swap(rawPoints[i][0], rawPoints[j][0]);
                swap(rawPoints[i][1], rawPoints[j][1]);
                swap(assocCurve[i], assocCurve[j]);
            }
        }
    }
    n_avg = 0;
    for (i = 0; i < n_pts; i++) {
        if (!visit[i]) {
            c = 0;
            px = 0;
            py = 0;
            for (j = 0; j < n_pts; j++) {
                if (!visit[j]) {
                    dx = rawPoints[i][0] - rawPoints[j][0];
                    dy = rawPoints[i][1] - rawPoints[j][1];
                    dm = dx * dx + dy * dy;
                    if (dm <= distsq &&
                        isInside(*assocCurve[i], rawPoints[j][0],
                                 rawPoints[j][1]) &&
                        isInside(*assocCurve[j], rawPoints[i][0],
                                 rawPoints[i][1])) {
                        visit[j] = true;
                        px += rawPoints[j][0];
                        py += rawPoints[j][1];
                        c++;
                    }
                }
            }
            if (c >= minPts) {
                avgPoints[n_avg][0] = px / c;
                avgPoints[n_avg][1] = py / c;
                n_avg++;
            }
        }
    }

    // Output results
    vector<cv::Point> pts(n_avg);
    n_out = 0;
    for (i = 0; i < n_avg; i++) {
        x = (int)avgPoints[i][0];
        y = (int)avgPoints[i][1];
        if (x > 0 && y > 0 && x < dp.ces_lo.ces[0].w &&
            y < dp.ces_lo.ces[0].h) {
            pts[n_out].x = x;
            pts[n_out].y = y;
            n_out++;
        }
    }

    // Free memory
    delete[] rawPoints;
    delete[] avgPoints;
    delete[] assocCurve;
    delete[] neighbors;

    return pts;
}

void *threadPhase2(void *t_data) {
    t_param *p = static_cast<t_param *>(t_data);
    int ns;
    ns = inflationSnake25d(*(p->dp), p->ix, p->iy, p->start_z, p->end_z,
                           SN25D_W_TENSION, SN25D_W_CURVATURE, SN25D_W_ZTENSION,
                           SN25D_W_ZCURVATURE, SN25D_W_ECURVE, SN25D_W_EINF_MIN,
                           SN25D_W_EINF_MAX, SN25D_W_EINF_STEP, p->ar_z,
                           *(p->outputSnakeList));

    pthread_mutex_lock(&tnslock);
    *(p->returnValue) += ns;
    pthread_mutex_unlock(&tnslock);
    sem_post(&sem2);
    pthread_exit(nullptr);
}

void retrieveAllSnakes25d(dataPacket &dp, vector<cv::Point> &initPts,
                          int start_z, int end_z,
                          snake25dList &outputSnakeList) {
    double ar_z = (int)LE_SSIZE_LO * TFACTOR / RESOLUTION;
    double ix, iy;
    volatile int tns = 0;
    uint n = initPts.size();

    int t;
    int rc;
    int numThreads = n;
    pthread_t *threads = new pthread_t[numThreads];
    pthread_attr_t attr;
    t_param *t_data = new t_param[numThreads];
    sem_init(&sem2, 0, numCores);
    pthread_mutex_unlock(&addlock);
    pthread_mutex_unlock(&tnslock);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE * 1024 * 1024);

    // Report
    cout << boost::format("Using %d initial points (slices: %04d - %04d)...") %
                n % start_z % end_z
         << endl;

    for (uint i = 0; i < n; i++) {
        sem_wait(&sem2);

        ix = (double)initPts[i].x;
        iy = (double)initPts[i].y;

        // Set thread params
        t_data[i].dp = &dp;
        t_data[i].ix = ix;
        t_data[i].iy = iy;
        t_data[i].start_z = start_z;
        t_data[i].end_z = end_z;
        t_data[i].ar_z = ar_z;
        t_data[i].outputSnakeList = &outputSnakeList;
        t_data[i].returnValue = &tns;
        //

        cout << "\rProcessing: " << (i + 1) << " / " << n
             << " (Total # of snakes: " << tns << ")";
        flush(cout);

        rc = pthread_create(&threads[i], &attr, threadPhase2,
                            static_cast<void *>(&t_data[i]));
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

    cout << endl << "Total # of snakes: " << tns << endl;

    delete[] threads;
    delete[] t_data;
}

void saveSnakeList(snake25dList &sarray, int tag_start_z, int tag_end_z) {
    string outputfn;
    outputfn = (boost::format("%s%s%d_%s.dat") % DESTPATH % "snakes25d_" %
                (tag_end_z - tag_start_z + 1) % FNAME)
                   .str();
    outputfn = (boost::format(outputfn) % tag_start_z).str();
    cout << "Saving: " << outputfn << endl;
    ofstream f(outputfn, ios::binary);

    uint n = sarray.size();
    f.write(reinterpret_cast<char *>(&n), sizeof(int));
    for (uint i = 0; i < n; i++) {
        f.write(reinterpret_cast<char *>(sarray[i].node),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d1_xy),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d2_xy),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d4_xy),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d1_z),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d2_z),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].d4_z),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].grad_Eint),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].grad_Eext),
                sizeof(double[SN_N][2]) * sn25d_t);
        f.write(reinterpret_cast<char *>(sarray[i].area),
                sizeof(double) * sn25d_t);
        f.write(reinterpret_cast<char *>(&sarray[i].status), sizeof(char));
        f.write(reinterpret_cast<char *>(&sarray[i].isValid),
                sizeof(char) * sn25d_t);
        f.write(reinterpret_cast<char *>(&sarray[i].validity), sizeof(double));
        f.write(reinterpret_cast<char *>(&sarray[i].start_z), sizeof(int));
        f.write(reinterpret_cast<char *>(&sarray[i].end_z), sizeof(int));
        f.write(reinterpret_cast<char *>(&sarray[i].t), sizeof(int));
        f.write(reinterpret_cast<char *>(&sarray[i].aspectRatio_z),
                sizeof(double));
    }
    f.close();
}

snake25dList loadSnakeArray(int tag_start_z, int tag_end_z) {
    int n;
    vector<snake25d> sarray;
    string inputfn;
    inputfn = (boost::format("%s%s%d_%s.dat") % DESTPATH % "snakes25d_" %
               (tag_end_z - tag_start_z + 1) % FNAME)
                  .str();
    inputfn = (boost::format(inputfn) % tag_start_z).str();
    cout << "Loading: " << inputfn << endl;
    ifstream f(inputfn, ios::binary);
    if (f.is_open()) {
        if (f.read(reinterpret_cast<char *>(&n), sizeof(int))) {
            if (n > 0) {
                sarray.resize(n);
                for (int i = 0; i < n; i++) {
                    f.read(reinterpret_cast<char *>(sarray[i].node),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d1_xy),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d2_xy),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d4_xy),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d1_z),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d2_z),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].d4_z),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].grad_Eint),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].grad_Eext),
                           sizeof(double[SN_N][2]) * sn25d_t);
                    f.read(reinterpret_cast<char *>(sarray[i].area),
                           sizeof(double) * sn25d_t);
                    f.read(reinterpret_cast<char *>(&(sarray[i].status)),
                           sizeof(char));
                    f.read(reinterpret_cast<char *>(sarray[i].isValid),
                           sizeof(char) * sn25d_t);
                    f.read(reinterpret_cast<char *>(&(sarray[i].validity)),
                           sizeof(double));
                    f.read(reinterpret_cast<char *>(&(sarray[i].start_z)),
                           sizeof(int));
                    f.read(reinterpret_cast<char *>(&(sarray[i].end_z)),
                           sizeof(int));
                    f.read(reinterpret_cast<char *>(&(sarray[i].t)),
                           sizeof(int));
                    f.read(reinterpret_cast<char *>(&(sarray[i].aspectRatio_z)),
                           sizeof(double));
                }
            }
        } else
            cerr << "Error loading data!" << endl;
    } else
        cerr << "Error loading data!" << endl;
    f.close();
    return sarray;
}

void savePoints(vector<cv::Point> &pts, int tag_start_z, int tag_end_z) {
    string outputfn;
    outputfn = (boost::format("%s%s%d_%s.dat") % DESTPATH % "initPts_" %
                (tag_end_z - tag_start_z + 1) % FNAME)
                   .str();
    outputfn = (boost::format(outputfn) % tag_start_z).str();
    cout << "Saving: " << outputfn << endl;
    ofstream f(outputfn, ios::binary);

    uint n = pts.size();
    f.write(reinterpret_cast<char *>(&n), sizeof(int));
    f.write(reinterpret_cast<char *>(pts.data()), sizeof(cv::Point) * n);

    f.close();
}

vector<cv::Point> loadPoints(int tag_start_z, int tag_end_z) {
    vector<cv::Point> pts;
    int n;
    string inputfn;
    inputfn = (boost::format("%s%s%d_%s.dat") % DESTPATH % "initPts_" %
               (tag_end_z - tag_start_z + 1) % FNAME)
                  .str();
    inputfn = (boost::format(inputfn) % tag_start_z).str();
    cout << "Loading: " << inputfn << endl;
    ifstream f(inputfn, ios::binary);
    if (f.is_open()) {
        f.read(reinterpret_cast<char *>(&n), sizeof(int));
        if (n > 0) {
            pts.resize(n);
            f.read(reinterpret_cast<char *>(pts.data()), sizeof(cv::Point) * n);
        }
        if (!f)
            cerr << "Error loading data!" << endl;
    } else
        cerr << "Error loading data!" << endl;
    f.close();
    return pts;
}

snake25dList filterSnakeListByValidity(snake25dList &in, double th) {
    snake25dList out;
    for (auto &s : in) {
        if (s.validity >= th) {
            out.push_back(s);
        }
    }
    return out;
}

vector<gradVec> getEnergyGradient(const vector<double> &ce, const int cw,
                                  const int ch) {
    int x, y;
    int xp, yp;
    double e_xp, e_yp, e;
    vector<gradVec> grad_ce(cw * ch);
    for (y = 0; y < ch; y++) {
        yp = y + 1;
        for (x = 0; x < cw; x++) {
            e = ce[y * cw + x];
            xp = x + 1;
            e_xp = e_yp = 0;
            if (xp < cw)
                e_xp = ce[y * cw + xp];
            if (yp < ch)
                e_yp = ce[yp * cw + x];
            grad_ce[y * cw + x].x = e_xp - e;
            grad_ce[y * cw + x].y = e_yp - e;
        }
    }
    return grad_ce;
}

vector<double> getCurveEnergyImage(vector<curve> &clist, const int w,
                                   const int h) {
    uint i;
    int j;
    curvePoints<double> points;
    cv::Mat ground;
    double max1, max2, min1, min2;
    int msize;
    int xi, yi;

    msize = w * h;

    // Alloc mem
    ground = cv::Mat::zeros(h, w, CV_32FC1);

    // Compute ground curve energy
    for (i = 0; i < clist.size(); i++) {
        points = getCurvePoints(clist[i], 1.0);

        for (j = 0; j < clist[i].len; j++) {
            xi = cvRound(points[j].x);
            yi = cvRound(points[j].y);
            if (xi >= 0 && yi >= 0 && xi < w && yi < h)
                ground.at<float>(yi, xi) += clist[i].avgscore;
        }
    }

    // Smoothing
    cv::minMaxLoc(ground, &min1, &max1);
    cv::GaussianBlur(ground, ground, cv::Size(0, 0), SN_GAUSSIAN, 0, 0);
    cv::minMaxLoc(ground, &min2, &max2);
    ground.convertTo(ground, -1, max1 / max2);

    // Alloc mem
    vector<double> ce(msize);

    // Copy to output array
    for (i = 0, yi = 0; yi < h; yi++) {
        for (xi = 0; xi < w; xi++, i++) {
            ce[i] = ground.at<float>(yi, xi);
        }
    }

    return ce;
}
