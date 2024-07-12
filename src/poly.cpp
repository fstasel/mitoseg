/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "poly.h"

void initPoly2d(double (*v)[2], int n, poly2d &p) {
    poly2d t(1);
    t[0].resize(n);

    for (int i = 0; i < n; i++) {
        t[0][i].x = (int)v[i][0];
        t[0][i].y = (int)v[i][1];
    }
    p.resize(1);
    combinePoly2d(t, poly2d(0), 0, p);
}

void initBlankPoly25d(int t, int start_z, poly25d &p) {
    p.t = t;
    p.start_z = start_z;
    p.end_z = start_z + t - 1;
    p.slice.clear();
    if (t)
        p.slice.resize(t);
}

void deinitPoly25d(poly25d &p) {
    p.slice.clear();
    p.t = 0;
    p.start_z = p.end_z = -1;
}

void initPoly25d(snake25d &s, poly25d &p) {
    initBlankPoly25d(s.t, s.start_z, p);
    for (int i = 0; i < s.t; i++) {
        initPoly2d(s.node[i], SN_N, p.slice[i]);
    }
}

vector<poly25d> initPoly25dArray(snake25dList &sarray) {
    vector<poly25d> parray(sarray.size());
    for (uint i = 0; i < sarray.size(); i++) {
        initPoly25d(sarray[i], parray[i]);
    }
    return parray;
}

void deinitPoly25dArray(vector<poly25d> &parray) {
    for (uint i = 0; i < parray.size(); i++) {
        deinitPoly25d(parray[i]);
    }
}

int poly2dArea(polygon &p, int *r) {
    int k, j;
    int area = 0;
    int n = p.size();
    for (k = 0, j = 1; k < n; k++, j = (j + 1) % n) {
        area += p[k].x * p[j].y - p[k].y * p[j].x;
    }
    area /= 2;
    if (area < 0) {
        if (r)
            *r = 1;
        area = -area;
    } else {
        if (r)
            *r = 0;
    }
    return area;
}

void correctPoly2d(poly2d &poly) {
    int i, j, n, r;

    for (auto &p : poly) {
        r = 0;
        poly2dArea(p, &r);
        if (r) {
            n = p.size();
            for (i = 0, j = n - 1; i < n / 2; i++, j--) {
                swap(p[i], p[j]);
            }
        }
    }
}

void combinePoly2d(const poly2d &p1, const poly2d &p2, int intersect,
                   poly2d &out, int *a1, int *a2, int *a12) {
    int w = 0, h = 0;
    int i;
    polygon k1;
    polygon k2;
    int n1, n2;

    poly2d::const_iterator p1_iter;
    poly2d::const_iterator p2_iter;

    // Determine working area
    p1_iter = p1.begin();
    p2_iter = p2.begin();
    do {
        if (p1_iter != p1.end() && p1_iter->size() > 0) {
            n1 = p1_iter->size();
            k1 = *p1_iter;
        } else {
            n1 = 0;
            k1.clear();
        }
        if (p2_iter != p2.end() && p2_iter->size() > 0) {
            n2 = p2_iter->size();
            k2 = *p2_iter;
        } else {
            n2 = 0;
            k2.clear();
        }

        for (i = 0; i < n1; i++) {
            if (w < k1[i].x)
                w = k1[i].x;
            if (h < k1[i].y)
                h = k1[i].y;
        }
        for (i = 0; i < n2; i++) {
            if (w < k2[i].x)
                w = k2[i].x;
            if (h < k2[i].y)
                h = k2[i].y;
        }

        // Go to next poly
        if (p1_iter != p1.end())
            p1_iter++;
        if (p2_iter != p2.end())
            p2_iter++;
    } while (p1_iter != p1.end() || p2_iter != p2.end());
    w += 2;
    h += 2;

    // Initialize buffers
    cv::Mat buf1 = cv::Mat::zeros(h, w, CV_8UC1);
    cv::Mat buf2 = cv::Mat::zeros(h, w, CV_8UC1);

    // Draw polygons
    p1_iter = p1.begin();
    p2_iter = p2.begin();
    do {
        if (p1_iter != p1.end() && p1_iter->size() > 0) {
            n1 = p1_iter->size();
            k1 = *p1_iter;
        } else {
            n1 = 0;
            k1.clear();
        }
        if (p2_iter != p2.end() && p2_iter->size() > 0) {
            n2 = p2_iter->size();
            k2 = *p2_iter;
        } else {
            n2 = 0;
            k2.clear();
        }

        if (n1)
            cv::fillPoly(buf1, (const cv::Point **)&k1, &n1, 1,
                         cv::Scalar(255, 255, 255));
        if (n2)
            cv::fillPoly(buf2, (const cv::Point **)&k2, &n2, 1,
                         cv::Scalar(255, 255, 255));

        // Go to next poly
        if (p1_iter != p1.end())
            p1_iter++;
        if (p2_iter != p2.end())
            p2_iter++;
    } while (p1_iter != p1.end() || p2_iter != p2.end());

    // Obtain area of the polygons
    if (a1)
        *a1 = cv::countNonZero(buf1);
    if (a2)
        *a2 = cv::countNonZero(buf2);

    // Create combined polygon
    if (intersect)
        cv::bitwise_and(buf1, buf2, buf1);
    else
        cv::bitwise_or(buf1, buf2, buf1);

    // Obtain combined polygon area
    if (a12)
        *a12 = cv::countNonZero(buf1);

    // Extract contours
    cv::findContours(buf1, out, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);
    correctPoly2d(out);
}

void combinePoly25d(const poly25d &p1, const poly25d &p2, int intersect,
                    poly25d &out, int *a1, int *a2, int *a12) {
    int start_z =
        (((p1.start_z < p2.start_z) & 1) ^ intersect) ? p1.start_z : p2.start_z;
    int end_z = (((p1.end_z > p2.end_z) & 1) ^ intersect) ? p1.end_z : p2.end_z;

    const poly2d *s1, *s2;
    poly2d temp;
    int t = end_z - start_z + 1;
    int z;
    int b1, b2, b12, c1 = 0, c2 = 0, c12 = 0;

    deinitPoly25d(out);
    if (t > 0) {
        initBlankPoly25d(t, start_z, out);

        for (z = start_z; z <= end_z; z++) {
            if (p1.start_z <= z && p1.end_z >= z)
                s1 = &p1.slice[z - p1.start_z];
            else
                s1 = &temp;
            if (p2.start_z <= z && p2.end_z >= z)
                s2 = &p2.slice[z - p2.start_z];
            else
                s2 = &temp;

            combinePoly2d(*s1, *s2, intersect, out.slice[z - start_z], &b1, &b2,
                          &b12);
            c1 += b1;
            c2 += b2;
            c12 += b12;
        }
    }

    if (a1)
        *a1 = c1;
    if (a2)
        *a2 = c2;
    if (a12)
        *a12 = c12;
}

vector<poly25d> mergeArrayOfPoly25d(vector<poly25d> &parray, double th) {
    uint i, j;
    uint n_out = 0;
    int a1, a2, a12;
    bool cont_merge = true;
    double d1, d2;

    poly25d temp;
    vector<poly25d> outArray;

    if (!parray.empty()) {
        outArray = parray;
        n_out = outArray.size();

        while (cont_merge) {
            cont_merge = false;

            for (i = 0; i < n_out - 1; i++) {
                cout << boost::format("\rMerging %d / %d (thres "
                                      "= %f)          ") %
                            i % (n_out - 1) % th;
                flush(cout);
                for (j = i + 1; j < n_out; j++) {
                    a1 = a2 = a12 = 0;
                    combinePoly25d(outArray[i], outArray[j], 1, temp, &a1, &a2,
                                   &a12);
                    if (a1)
                        d1 = (double)a12 / a1;
                    else
                        d1 = 0;
                    if (a2)
                        d2 = (double)a12 / a2;
                    else
                        d2 = 0;
                    // Merging condition
                    if (d1 >= th || d2 >= th) {
                        combinePoly25d(outArray[i], outArray[j], 0, temp);
                        outArray[j] = outArray[n_out - 1];
                        outArray[i] = temp;
                        deinitPoly25d(outArray[n_out - 1]);
                        n_out--;
                        cont_merge = true;
                    }
                }
            }
        }
        cout << endl;
    }

    outArray.resize(n_out);
    return outArray;
}

vector<poly25d> convertValidSnakesToPolyArray(snake25dList &sarray, double th) {
    snake25dList vs = filterSnakeListByValidity(sarray, th);
    return initPoly25dArray(vs);
}

int getTriangles(polygon &p1, polygon &p2, vector<int3> &faces) {
    faces.clear();
    faces.resize(p1.size() + p2.size());
    uint n_faces = 0;

    uint i, j, ii, jj, k;
    int dx, dy, d;
    int dmin = INT_MAX;
    uint is = 0, js = 0;

    for (i = 0; i < p1.size(); i++) {
        for (j = 0; j < p2.size(); j++) {
            d = 0;
            ii = i;
            jj = j;
            for (k = 0; k < 20; k++) {
                dx = p1[ii].x - p2[jj].x;
                dy = p1[ii].y - p2[jj].y;
                d += 1 + dx * dx + dy * dy;
                ii++;
                if (ii >= p1.size())
                    ii -= p1.size();
                jj++;
                if (jj >= p2.size())
                    jj -= p2.size();
            }
            if (d < dmin) {
                dmin = d;
                is = i + 10;
                if (is >= p1.size())
                    is -= p1.size();
                js = j + 10;
                if (js >= p2.size())
                    js -= p2.size();
            }
        }
    }

    double id, jd;
    i = j = 0;
    while (n_faces < p1.size() + p2.size()) {
        id = (double)i / p1.size();
        jd = (double)j / p2.size();
        ii = (i + is) % p1.size();
        jj = (j + js) % p2.size();
        if (id < jd) {
            faces[n_faces].index[0] = ii;
            faces[n_faces].index[1] = (ii + 1) % p1.size();
            faces[n_faces].index[2] = jj + p1.size();
            i++;
        } else {
            faces[n_faces].index[0] = ii;
            faces[n_faces].index[2] = jj + p1.size();
            faces[n_faces].index[1] = (jj + 1) % p2.size() + p1.size();
            j++;
        }
        n_faces++;
    }

    return n_faces;
}

void savePolyArrayAsPLY(vector<poly25d> &p) {
    double xy_scale = LE_SSIZE_LO * TFACTOR;
    double z_scale = RESOLUTION;
    double x, y, z;
    uint n = p.size();
    uint i, j, k, n_vertex = 0, cv;
    ofstream f;
    vector<int3> faces;
    uint n_faces;

    for (k = 0; k < n; k++) {
        for (i = 0; i < p[k].t; i++) {
            for (auto &s : p[k].slice[i])
                n_vertex += s.size();
        }
    }

    n_faces = 0;
    for (k = 0; k < n; k++) {
        for (i = 0; i < p[k].t; i++) {
            if (i < p[k].t - 1) {
                n_faces +=
                    getTriangles(p[k].slice[i][0], p[k].slice[i + 1][0], faces);
            }
        }
    }

    string outputfn;
    outputfn = (boost::format("%s%s%s.ply") % DESTPATH % "poly_" % FNAME).str();
    outputfn = (boost::format(outputfn) % SLICE_START).str();
    cout << "Saving: " << outputfn << endl;
    f.open(outputfn);

    f << "ply" << endl
      << "format ascii 1.0" << endl
      << "element vertex " << n_vertex << endl
      << "property float x" << endl
      << "property float y" << endl
      << "property float z" << endl
      << "element face " << n_faces << endl
      << "property list uchar int vertex_indices" << endl
      << "end_header" << endl;

    for (k = 0; k < n; k++) {
        for (i = 0; i < p[k].t; i++) {
            z = (i + p[k].start_z) * z_scale;
            for (auto &s : p[k].slice[i]) {
                for (j = 0; j < s.size(); j++) {
                    x = s[j].x * xy_scale;
                    y = s[j].y * xy_scale;
                    f << boost::format("%f %f %f") % x % y % z << endl;
                }
            }
        }
    }

    cv = 0;
    for (k = 0; k < n; k++) {
        for (i = 0; i < p[k].t; i++) {
            if (i < p[k].t - 1) {
                n_faces =
                    getTriangles(p[k].slice[i][0], p[k].slice[i + 1][0], faces);
                for (j = 0; j < n_faces; j++) {
                    f << boost::format("3 %d %d %d") %
                             (faces[j].index[0] + cv) %
                             (faces[j].index[1] + cv) % (faces[j].index[2] + cv)
                      << endl;
                }
            }
            for (auto &s : p[k].slice[i]) {
                cv += s.size();
                if (p[k].slice[i].size() > 1)
                    cout << "Warning: Poly2d list is not supported!" << endl;
            }
        }
    }

    f.close();
}
