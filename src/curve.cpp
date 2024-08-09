/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "curve.h"

#include <boost/format.hpp>

vector<double> angle_x(LE_BINS);
vector<double> angle_y(LE_BINS);

vector<localEnergy> getLocalEnergyMap(vector<ridge> &r, const int width,
                                      const int height, const int ssize,
                                      const int bsize, int *w, int *h) {
    *w = (width / ssize) + ((width % ssize) ? 1 : 0);
    *h = (height / ssize) + ((height % ssize) ? 1 : 0);
    int maxi = (*w) * (*h);
    int maxj = width * height;
    int hsize = bsize / 2;
    int x, y, x1, y1, x2, y2, xx, yy, i, j, b, major;
    float m;
    vector<localEnergy> e(maxi);

    // Init bin indices
    int *q = new int[maxj];
    for (j = 0; j < maxj; j++) {
        q[j] = cvRound(LE_BINS * r[j].angle / 180.0f) % LE_BINS;
    }

    // Main loop
    for (y = 0, i = 0; y < height; y += ssize) {
        y1 = y - hsize;
        y2 = y + hsize;
        if (y1 < 0)
            y1 = 0;
        if (y2 > height)
            y2 = height;

        for (x = 0; x < width; x += ssize, i++) {
            x1 = x - hsize;
            x2 = x + hsize;
            if (x1 < 0)
                x1 = 0;
            if (x2 > width)
                x2 = width;

            // Calculate energy per direction
            for (yy = y1; yy < y2; yy++) {
                for (xx = x1; xx < x2; xx++) {
                    j = yy * width + xx;
                    if (r[j].mask) {
                        e[i].energy[q[j]] += r[j].mag;
                    }
                }
            }

            // Find major direction
            major = 0;
            m = 0;
            for (b = 0; b < LE_BINS; b++) {
                if (e[i].energy[b] > m) {
                    m = e[i].energy[b];
                    major = b;
                }
            }
            e[i].major = major;
        }
    }

    // Free mem.
    delete[] q;

    return e;
}

cv::Mat visualizeLEMajority(vector<localEnergy> &e, const int w, const int h,
                            const int ssize, const int showDir) {
    float maxe = 0;
    double a;
    int i, x, y, b, k, xx, yy, mx, my, x1, x2, y1, y2, dls;

    cv::Mat visual(ssize * h, ssize * w, CV_8UC1);

    dls = ssize / 2 - 2;
    if (dls < 1)
        dls = 1;

    // Find global highest energy
    for (y = 0, i = 0; y < h; y++) {
        for (x = 0; x < w; x++, i++) {
            for (b = 0; b < LE_BINS; b++) {
                if (e[i].energy[b] > maxe)
                    maxe = e[i].energy[b];
            }
        }
    }

    // Main loop
    for (y = 0, i = 0; y < h; y++) {
        yy = y * ssize;

        for (x = 0; x < w; x++, i++) {
            xx = x * ssize;
            k = cvRound(255 * e[i].energy[e[i].major] / maxe);
            cv::rectangle(visual, cv::Point(xx, yy),
                          cv::Point(xx + ssize - 1, yy + ssize - 1),
                          cv::Scalar(k), cv::FILLED);
            if (showDir) {
                a = (PI / 2) - (PI * e[i].major) / LE_BINS;
                mx = cvRound(cos(a) * dls);
                my = cvRound(sin(a) * dls);
                x1 = x2 = xx + ssize / 2;
                y1 = y2 = yy + ssize / 2;
                x1 += mx;
                y1 += my;
                x2 -= mx;
                y2 -= my;
                if (k > 128)
                    k = 0;
                else
                    k = 255;
                cv::line(visual, cv::Point(x1, y1), cv::Point(x2, y2),
                         cv::Scalar(k));
            }
        }
    }
    return visual;
}

void calcCurveParams(vector<localEnergy> &e, const int w, const int h, curve &c,
                     const double S) {
    double t, k, dt, m, wtheta;
    double x, y, xp, yp, xn, yn, mx, my;
    int xi, yi, i, j;

    // Calculate model params
    int dx = c.x2 - c.x1;
    int dy = c.y2 - c.y1;
    double R = sqrt((double)(dx * dx + dy * dy));
    double b = (double)(4 * c.h) / R;
    double a = -b / R;
    double cost = dx / R;
    double sint = dy / R;

    c.R = R;
    c.a = a;
    c.b = b;
    c.sint = sint;
    c.cost = cost;
    c.score = 0;
    c.len = 0;
    if (R < EPS)
        return;

    // Main loop
    t = -S / sqrt(1.0 + b * b);
    k = a * t * t + b * t;
    xp = cost * t - sint * k + c.x1;
    yp = sint * t + cost * k + c.y1;
    x = (double)c.x1;
    y = (double)c.y1;
    dt = -t;
    t = 0;
    while (true) {
        t += dt;
        k = a * t * t + b * t;

        xn = cost * t - sint * k + c.x1;
        yn = sint * t + cost * k + c.y1;

        mx = xn - xp;
        my = yn - yp;

        xi = cvRound(x);
        yi = cvRound(y);
        if (xi >= 0 && yi >= 0 && xi < w && yi < h) {
            i = yi * w + xi;
            for (j = 0; j < LE_BINS; j++) {
                wtheta = angle_x[j] * my + angle_y[j] * mx;
                wtheta = 2 * (wtheta * wtheta) / (mx * mx + my * my) - 1;
                c.score += (float)wtheta * e[i].energy[j];
            }
        }

        c.len++;

        if (t > R)
            break;

        xp = x;
        yp = y;
        x = xn;
        y = yn;
        m = 2 * a * t + b;
        dt = S / sqrt(1.0 + m * m);
    }
}

void drawBinaryCurve(cv::Mat &im, curve &c, const int ssize, const int psize,
                     const double S, cv::Scalar cl) {
    double t, k, dt, m, x, y;
    int hsize = ssize / 2;

    t = 0;
    while (t <= c.R && c.R >= EPS) {
        k = c.a * t * t + c.b * t;
        x = c.cost * t - c.sint * k + c.x1;
        y = c.sint * t + c.cost * k + c.y1;

        cv::circle(
            im,
            cv::Point(cvRound(x * ssize) + hsize, cvRound(y * ssize) + hsize),
            psize, cl, cv::FILLED);

        m = 2 * c.a * t + c.b;
        dt = S / sqrt(1.0 + m * m);
        t += dt;
    }
}

void drawCurve(cv::Mat &im, curve &c, const int ssize, const double S,
               cv::Scalar cl) {
    double t, k, dt, m, x, y;
    int hsize = 0;
    cv::Point p1, p2;
    bool first = true;

    t = 0;
    while (t <= c.R && c.R >= EPS) {
        k = c.a * t * t + c.b * t;
        x = c.cost * t - c.sint * k + c.x1;
        y = c.sint * t + c.cost * k + c.y1;

        if (first) {
            first = false;
            p2 = cv::Point(cvRound(x * ssize) + hsize,
                           cvRound(y * ssize) + hsize);
        } else {
            p1 = p2;
            p2 = cv::Point(cvRound(x * ssize) + hsize,
                           cvRound(y * ssize) + hsize);
            cv::line(im, p1, p2, cl, 1);
        }

        m = 2 * c.a * t + c.b;
        dt = S / sqrt(1.0 + m * m);
        t += dt;
    }
}

int fitCurve(vector<localEnergy> &e, const int w, const int h, curve &c,
             const int maxIter, char *mask) {
    int i, dx, dy;
    curve best, cur, max, n;

    const double maxcurv = FC_MAXCURV * FC_MAXCURV;

    // Initial conditions
    int proc = 0;
    int iter = 0;
    int fail = 1;

    // Calculate score of initial curve
    cur = c;
    calcCurveParams(e, w, h, cur, FC_STEP);
    best = cur;

    while (proc < 2 && iter < maxIter) {
        // Use max as current
        max = cur;
        n = cur;

        if (proc == 0) {
            // Find best neighbor for (x1, y1) for all neighboring y1
            for (n.y1 = cur.y1 - FC_XYRANGE; n.y1 <= cur.y1 + FC_XYRANGE;
                 n.y1 += FC_XYSTEP) {
                // if inside range
                if (n.y1 >= 0 && n.y1 < h) {
                    dy = n.y2 - n.y1;
                    dy = dy * dy;
                    // for all neighboring x1
                    for (n.x1 = cur.x1 - FC_XYRANGE;
                         n.x1 <= cur.x1 + FC_XYRANGE; n.x1 += FC_XYSTEP) {
                        // if inside range
                        if (n.x1 >= 0 && n.x1 < w) {
                            dx = n.x2 - n.x1;
                            dx = dx * dx;

                            // If not masked
                            i = n.y1 * w + n.x1;
                            if (!(mask && mask[i])) {
                                // for all neighboring h
                                for (n.h = cur.h - FC_HRANGE;
                                     n.h <= cur.h + FC_HRANGE;
                                     n.h += FC_HSTEP) {
                                    // If curvature condition is satisfied
                                    if (n.h * n.h < maxcurv * (dx + dy)) {
                                        // Calculate score
                                        calcCurveParams(e, w, h, n, FC_STEP);
                                        if (n.score > max.score) {
                                            max = n;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else if (proc == 1) {
            // Find best neighbor for (x2, y2) for all neighboring y2
            for (n.y2 = cur.y2 - FC_XYRANGE; n.y2 <= cur.y2 + FC_XYRANGE;
                 n.y2 += FC_XYSTEP) {
                // if inside range
                if (n.y2 >= 0 && n.y2 < h) {
                    dy = n.y2 - n.y1;
                    dy = dy * dy;

                    // for all neighboring x2
                    for (n.x2 = cur.x2 - FC_XYRANGE;
                         n.x2 <= cur.x2 + FC_XYRANGE; n.x2 += FC_XYSTEP) {
                        // if inside range
                        if (n.x2 >= 0 && n.x2 < w) {
                            dx = n.x2 - n.x1;
                            dx = dx * dx;

                            // If not masked
                            i = n.y1 * w + n.x1;
                            if (!(mask && mask[i])) {
                                // for all neighboring h
                                for (n.h = cur.h - FC_HRANGE;
                                     n.h <= cur.h + FC_HRANGE;
                                     n.h += FC_HSTEP) {
                                    // If curvature condition is satisfied
                                    if (n.h * n.h < maxcurv * (dx + dy)) {
                                        // Calculate score
                                        calcCurveParams(e, w, h, n, FC_STEP);
                                        if (n.score > max.score) {
                                            max = n;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        cur = max;

        // Take current curve if better
        if (best.score < cur.score) {
            best = cur;
            fail = 0;
        } else {
            // Go on next step
            proc++;
        }

        // Increment iteration
        iter++;
    }

    if (!fail) {
        c = best;
    }

    return fail;
}

float energystr(localEnergy &e) {
    const int r = (LE_BINS / 4) - 1;
    float se = 0;
    for (int j = e.major - r; j <= e.major + r; j++) {
        se += e.energy[j];
    }
    return se;
}

curveList fitCurves(vector<localEnergy> &e, const int w, const int h,
                    const int FCS_XYSTEP, const int FCS_XYRANGE,
                    const float FCS_INIT_THRESH, const int FCS_MINLEN) {
    int nc = 0;
    int x, y, i, j, xi, yi, xm, ym, x1, x2, y1, y2;
    float sm, si, gm;
    curve c;
    curveList clist;
    bool *scanned = new bool[w * h]();

    // Find global max
    gm = 0;
    for (i = 0; i < w * h; i++) {
        si = energystr(e[i]);
        if (si > gm)
            gm = si;
    }

    // Scan blocks for local max
    for (y = 0; y < h; y += FCS_XYSTEP) {
        y1 = y + FCS_XYSTEP / 2 - FCS_XYRANGE;
        if (y1 < 0)
            y1 = 0;
        y2 = y + FCS_XYSTEP / 2 + FCS_XYRANGE;
        if (y2 > h)
            y2 = h;

        for (x = 0; x < w; x += FCS_XYSTEP) {
            x1 = x + FCS_XYSTEP / 2 - FCS_XYRANGE;
            if (x1 < 0)
                x1 = 0;
            x2 = x + FCS_XYSTEP / 2 + FCS_XYRANGE;
            if (x2 > w)
                x2 = w;

            xm = x1;
            ym = y1;
            i = ym * w + xm;
            sm = energystr(e[i]);
            for (yi = y1; yi < y2; yi++) {
                for (xi = x1; xi < x2; xi++) {
                    i = yi * w + xi;
                    si = energystr(e[i]);
                    if (si > sm) {
                        sm = si;
                        xm = xi;
                        ym = yi;
                    }
                }
            }

            i = ym * w + xm;
            // If (not visited) and strong enough
            if (!scanned[i] && sm / gm > FCS_INIT_THRESH) {
                scanned[i] = true;
                c.x1 = c.x2 = xm;
                c.y1 = c.y2 = ym;
                c.h = 0;
                if (!fitCurve(e, w, h, c, FC_MAXITER)) {
                    // If curve is long enough
                    if (c.len >= FCS_MINLEN) {
                        clist.push_back(c);
                        nc++;
                    }
                }
            }
        }
    }

    // Clean up curve array for duplicates
    for (i = 0; i < nc - 1; i++) {
        for (j = i + 1; j < nc; j++) {
            if (clist[i].x1 == clist[j].x1 && clist[i].x2 == clist[j].x2 &&
                clist[i].y1 == clist[j].y1 && clist[i].y2 == clist[j].y2 &&
                clist[i].h == clist[j].h) {
                clist[j] = clist[nc - 1];
                nc--;
                j--;
            }
        }
    }
    clist.resize(nc);

    // The other stuff
    for (i = 0; i < nc; i++) {
        // Correct sign of h param (make positive)
        if (clist[i].h < 0) {
            clist[i].h *= -1;
            clist[i].a *= -1;
            clist[i].b *= -1;
            clist[i].cost *= -1;
            clist[i].sint *= -1;
            x = clist[i].x1;
            y = clist[i].y1;
            clist[i].x1 = clist[i].x2;
            clist[i].y1 = clist[i].y2;
            clist[i].x2 = x;
            clist[i].y2 = y;
        }

        // Calculate aux. curve params
        clist[i].avgscore = clist[i].score / clist[i].len;
    }

    delete[] scanned;
    return clist;
}

void visualizeCurves(cv::Mat &bckgnd, curveList &clist, const int ssize,
                     cv::Mat &visual, const float FLC_THRESH,
                     const float FLC_AVGTHRESH) {
    cv::Scalar cl;
    vector<cv::Mat> channels(3, bckgnd);
    cv::merge(channels, visual);

    for (uint i = 0; i < clist.size(); i++) {
        if (clist[i].score > FLC_THRESH && clist[i].avgscore > FLC_AVGTHRESH) {
            cl = cv::Scalar(255, 0, 0);
        } else {
            cl = cv::Scalar(0, 0, 255);
        }
        drawCurve(visual, clist[i], ssize, FC_STEP, cl);
    }
}

void initCurveFitting() {
    int j;
    double t;
    for (j = 0; j < LE_BINS; j++) {
        t = (PI * j) / LE_BINS;
        angle_x[j] = cos(t);
        angle_y[j] = sin(t);
    }
}

cv::Mat drawBinaryCurves(curveList &clist, const int psize, const int w,
                         const int h) {
    int r = psize / 2;
    cv::Mat out = cv::Mat::zeros(h, w, CV_8UC1);

    for (uint i = 0; i < clist.size(); i++) {
        drawBinaryCurve(out, clist[i], 1, r, FC_STEP, cv::Scalar(255));
    }

    return out;
}

void visualizeHiLoCurves(curveList &clist_lo, const int w_lo, const int h_lo,
                         curveList &clist_hi, const double maxoverlap,
                         cv::Mat &bckgnd, cv::Mat &visual) {
    // Prepare low freq. curves coverage
    cv::Mat memtemp =
        drawBinaryCurves(clist_lo, (int)CURV_HILO_COV, w_lo, h_lo);

    vector<cv::Mat> channels(3, bckgnd);
    cv::merge(channels, visual);

    double scale = (double)LE_SSIZE_LO / LE_SSIZE_HI;
    int xi, yi, in, tot;
    uint i;
    double t, k, x, y, m, dt;

    for (i = 0; i < clist_hi.size(); i++) {
        curve &c = clist_hi[i];

        in = 0;
        tot = 0;

        t = 0;
        while (t <= c.R && c.R >= EPS) {
            k = c.a * t * t + c.b * t;
            x = c.cost * t - c.sint * k + c.x1;
            y = c.sint * t + c.cost * k + c.y1;

            xi = cvRound(x * scale);
            yi = cvRound(y * scale);
            if (xi >= 0 && yi >= 0 && xi < w_lo && yi < h_lo) {
                if (memtemp.at<unsigned char>(yi, xi))
                    in++;
            }
            tot++;

            m = 2 * c.a * t + c.b;
            dt = FC_STEP / sqrt(1.0 + m * m);
            t += dt;
        }

        if ((double)in / tot < maxoverlap) {
            drawCurve(visual, c, (int)LE_SSIZE_HI, FC_STEP,
                      cv::Scalar(0, 0, 255));
        }
    }

    for (i = 0; i < clist_lo.size(); i++) {
        drawCurve(visual, clist_lo[i], (int)LE_SSIZE_LO, FC_STEP,
                  cv::Scalar(255, 0, 0));
    }
}

curveList getHixCurves(curveList &clist_lo, const int w_lo, const int h_lo,
                       curveList &clist_hi) {
    curveList clist_hix;

    // Prepare low freq. curves coverage
    cv::Mat memtemp =
        drawBinaryCurves(clist_lo, (int)CURV_HILO_COV, w_lo, h_lo);

    double scale = (double)LE_SSIZE_LO / LE_SSIZE_HI;
    int xi, yi, in, tot;
    curve *c;
    double t, k, x, y, m, dt;

    for (uint i = 0; i < clist_hi.size(); i++) {
        c = &clist_hi[i];

        in = 0;
        tot = 0;

        t = 0;
        while (t <= c->R && c->R >= EPS) {
            k = c->a * t * t + c->b * t;
            x = c->cost * t - c->sint * k + c->x1;
            y = c->sint * t + c->cost * k + c->y1;

            xi = cvRound(x * scale);
            yi = cvRound(y * scale);
            if (xi >= 0 && yi >= 0 && xi < w_lo && yi < h_lo) {
                if (memtemp.at<unsigned char>(yi, xi))
                    in++;
            }
            tot++;

            m = 2 * c->a * t + c->b;
            dt = FC_STEP / sqrt(1.0 + m * m);
            t += dt;
        }

        if ((double)in / tot < CURV_HILO_COV_TH) {
            clist_hix.push_back(clist_hi[i]);
        }
    }

    return clist_hix;
}

curveList filterCurves(const curveList &c_in, const float FLC_THRESH,
                       const float FLC_AVGTHRESH) {
    curveList c_out;
    for (uint i = 0; i < c_in.size(); i++) {
        if (c_in[i].score > FLC_THRESH && c_in[i].avgscore > FLC_AVGTHRESH) {
            c_out.push_back(c_in[i]);
        }
    }
    return c_out;
}

curvePoints<double> getCurvePoints(const curve &c, const double S) {
    int j;
    double t, k, dt, m, x, y;
    curvePoints<double> points;

    // Init point counter
    j = 0;
    // Scan "len" points (0<=t<=R)
    t = 0;
    while (j < c.len) {
        k = c.a * t * t + c.b * t;
        x = c.cost * t - c.sint * k + c.x1;
        y = c.sint * t + c.cost * k + c.y1;

        points.push_back({x, y});

        m = 2 * c.a * t + c.b;
        dt = S / sqrt(1.0 + m * m);
        t += dt;
        j++;
    }

    return points;
}

void saveLocalEnergyStructure(int s, const string tag, vector<localEnergy> &e,
                              int w, int h) {
    if (!e.empty()) {
        string outputfn;
        outputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
        outputfn = (boost::format(outputfn) % s).str();
        cout << "Saving: " << outputfn << endl;
        ofstream f(outputfn, ios::binary);
        f.write(reinterpret_cast<char *>(&w), sizeof(int));
        f.write(reinterpret_cast<char *>(&h), sizeof(int));
        f.write(reinterpret_cast<char *>(e.data()),
                sizeof(localEnergy) * w * h);
        f.close();
    }
}

void loadLocalEnergyStructure(int s, const string tag, vector<localEnergy> &e,
                              int *w, int *h) {
    string inputfn;
    inputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
    inputfn = (boost::format(inputfn) % s).str();
    e.clear();
    cout << "Loading: " << inputfn << endl;
    ifstream f(inputfn, ios::binary);
    f.read(reinterpret_cast<char *>(w), sizeof(int));
    f.read(reinterpret_cast<char *>(h), sizeof(int));
    e.resize((*w) * (*h));
    if (f.read(reinterpret_cast<char *>(e.data()),
               sizeof(localEnergy) * (*w) * (*h)))
        cerr << "Error loading data!" << endl;
    f.close();
}

void saveCurveStructure(int s, const string tag, curveList &c) {
    if (!c.empty()) {
        string outputfn;
        outputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
        outputfn = (boost::format(outputfn) % s).str();
        cout << "Saving: " << outputfn << endl;
        ofstream f(outputfn, ios::binary);
        uint numRecords = c.size();
        f.write(reinterpret_cast<char *>(&numRecords), sizeof(uint));
        f.write(reinterpret_cast<char *>(c.data()), sizeof(curve) * numRecords);
        f.close();
    }
}

curveList loadCurveStructure(int s, const string tag) {
    int numRecords;
    string inputfn;
    inputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
    inputfn = (boost::format(inputfn) % s).str();
    cout << "Loading: " << inputfn << endl;
    ifstream f(inputfn, ios::binary);
    f.read(reinterpret_cast<char *>(&numRecords), sizeof(int));
    curveList c(numRecords);
    if (!f.read(reinterpret_cast<char *>(c.data()), sizeof(curve) * numRecords))
        cerr << "Error loading data!" << endl;
    f.close();
    return c;
}

int isInside(curve &c, double x, double y, double bias) {
    x -= c.x1;
    y -= c.y1;
    double t = c.cost * x + c.sint * y;
    double k = -c.sint * x + c.cost * y;
    if (k + bias <= c.a * t * t + c.b * t)
        return 1;
    return 0;
}
