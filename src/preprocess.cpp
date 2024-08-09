/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "preprocess.h"

void autoLevels(cv::Mat &input, cv::Mat &output) {
    int i, a, b;
    double prob[256], p, c, m;
    int histSize[] = {256};
    float range[] = {0, 256};
    const float *ranges[] = {range};
    cv::Mat hist;
    cv::calcHist(&input, 1, 0, cv::Mat(), hist, 1, histSize, ranges);

    for (i = 0; i < 256; i++)
        prob[i] = (double)hist.at<float>(i) / (input.cols * input.rows);
    a = 0;
    p = 0;
    while (p < AUTOLEVELSCUT) {
        p += prob[a];
        a++;
    }
    b = 255;
    p = 0;
    while (p < AUTOLEVELSCUT) {
        p += prob[b];
        b--;
    }
    m = 255.0 / (b - a);
    c = -m * a;
    for (int y = 0; y < input.rows; y++) {
        for (int x = 0; x < input.cols; x++) {
            int i = cvRound(input.at<uchar>(y, x) * m + c);
            i = clamp(i, 0, 255);
            output.at<uchar>(y, x) = i;
        }
    }
}

cv::Mat getResampledImage(const cv::Mat &input) {
    double sf;
    cv::Mat output, temp;
    if (TFACTOR == RESOLUTION) {
        output = cv::Mat(input.cols, input.rows, CV_32FC1);
        input.convertTo(output, CV_32FC1);
    } else {
        sf = RESOLUTION / TFACTOR;
        input.convertTo(temp, CV_32FC1);
        cv::resize(temp, output,
                   cv::Size(cvRound(input.cols * sf), cvRound(input.rows * sf)),
                   0, 0, cv::INTER_LINEAR);
    }
    return output;
}

void getRidges(cv::Mat &input, vector<ridge> &r, cv::Mat &output,
               cv::Mat &mask) {
    cv::Mat image;
    cv::Mat gxx, gxy, gyx, gyy;
    cv::Mat V, E;

    int i, j, ind;
    float m[4], p, a;
    cv::Mat M;

    gxx.create(input.size(), CV_32FC1);
    gxy.create(input.size(), CV_32FC1);
    gyx.create(input.size(), CV_32FC1);
    gyy.create(input.size(), CV_32FC1);
    V.create(2, 2, CV_32F);
    E.create(2, 1, CV_32F);

    input.copyTo(image);

    // Ridge detection using Hessian
    cv::Sobel(image, gxx, CV_32FC1, 2, 0, 3, 1, 0, cv::BORDER_REPLICATE);
    cv::Sobel(image, gxy, CV_32FC1, 1, 1, 3, 1, 0, cv::BORDER_REPLICATE);
    cv::Sobel(image, gyx, CV_32FC1, 1, 1, 3, 1, 0, cv::BORDER_REPLICATE);
    cv::Sobel(image, gyy, CV_32FC1, 0, 2, 3, 1, 0, cv::BORDER_REPLICATE);

    for (i = 0, ind = 0; i < image.rows; i++) {
        for (j = 0; j < image.cols; j++, ind++) {
            m[0] = gxx.at<float>(i, j);
            m[1] = gxy.at<float>(i, j);
            m[2] = gyx.at<float>(i, j);
            m[3] = gyy.at<float>(i, j);

            M = cv::Mat(2, 2, CV_32F, m);
            cv::eigen(M, E, V);

            if (E.at<float>(0) < 0) {
                p = 0;
                mask.at<uchar>(i, j) = 0;
                r[ind].mask = false;
            } else {
                if (E.at<float>(1) < 0)
                    p = E.at<float>(0);
                else
                    p = E.at<float>(0) - E.at<float>(1);
                mask.at<uchar>(i, j) = 255;
                r[ind].mask = true;
            }

            r[ind].mag = p;

            a = (float)(atan2(-V.at<float>(0, 1), V.at<float>(0, 0)) * 180 /
                        PI);
            if (a < 0)
                a += 360;
            r[ind].angle = a;
            r[ind].vec[0] = V.at<float>(0, 0);
            r[ind].vec[1] = V.at<float>(0, 1);
            r[ind].val[0] = E.at<float>(0);
            r[ind].val[1] = E.at<float>(1);
            r[ind].dir = (a > 180) ? a - 180 : a;

            // dir4
            if ((a > 22.5 && a <= 67.5) || (a > 202.5 && a <= 247.5))
                r[ind].dir4 = 1;
            else if ((a > 67.5 && a <= 112.5) || (a > 247.5 && a <= 292.5))
                r[ind].dir4 = 2;
            else if ((a > 112.5 && a <= 157.5) || (a > 292.5 && a <= 337.5))
                r[ind].dir4 = 3;
            else
                r[ind].dir4 = 0;

            // dir8
            if ((a > 11.25 && a <= 33.75) || (a > 191.25 && a <= 213.75))
                r[ind].dir8 = 1;
            else if ((a > 33.75 && a <= 56.25) || (a > 213.75 && a <= 236.25))
                r[ind].dir8 = 2;
            else if ((a > 56.25 && a <= 78.75) || (a > 236.25 && a <= 258.75))
                r[ind].dir8 = 3;
            else if ((a > 78.75 && a <= 101.25) || (a > 258.75 && a <= 281.25))
                r[ind].dir8 = 4;
            else if ((a > 101.25 && a <= 123.75) || (a > 281.25 && a <= 303.75))
                r[ind].dir8 = 5;
            else if ((a > 123.75 && a <= 146.25) || (a > 303.75 && a <= 326.25))
                r[ind].dir8 = 6;
            else if ((a > 146.25 && a <= 168.75) || (a > 326.25 && a <= 348.75))
                r[ind].dir8 = 7;
            else
                r[ind].dir8 = 0;

            output.at<float>(i, j) = r[ind].mag;
        }
    }
}

void visualizeRidgesDir(int opt, vector<ridge> &r, cv::Mat &visual_dir) {
    int i, j, ind, mind = visual_dir.cols * visual_dir.rows, p;
    float ridges_max = r[0].mag, ridges_min = r[0].mag;

    for (ind = 0; ind < mind; ind++) {
        if (ridges_max < r[ind].mag)
            ridges_max = r[ind].mag;
        if (ridges_min > r[ind].mag)
            ridges_min = r[ind].mag;
    }

    for (i = 0, ind = 0; i < visual_dir.rows; i++) {
        for (j = 0; j < visual_dir.cols; j++, ind++) {
            p = cvRound(255.0f * (r[ind].mag - ridges_min) /
                        (ridges_max - ridges_min));

            if (opt == 8) {
                if (r[ind].dir8 == 7)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, 0, 0);
                else if (r[ind].dir8 == 6)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(0, p, 0);
                else if (r[ind].dir8 == 5)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, p);
                else if (r[ind].dir8 == 4)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, p, 0);
                else if (r[ind].dir8 == 3)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, 0, p);
                else if (r[ind].dir8 == 2)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(0, p, p);
                else if (r[ind].dir8 == 1)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, p, p);
                else
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, p, p / 2);
            } else if (opt == 4) {
                if (r[ind].dir4 == 3)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, 0, 0);
                else if (r[ind].dir4 == 2)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(0, p, 0);
                else if (r[ind].dir4 == 1)
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, p);
                else
                    visual_dir.at<cv::Vec3b>(i, j) = cv::Vec3b(p, p, p);
            }
        }
    }
}

void normalizeRidges(vector<ridge> &r, cv::Mat &ridges, cv::Mat &mask) {
    // 0.5% cut for norm.dist.
    const double autoLevelCut = 2.81;

    cv::Scalar mean, std;
    float p;
    int x, y, i;
    cv::meanStdDev(ridges, mean, std, mask);
    double a, b, m, c;

    a = mean.val[0] - std.val[0] * autoLevelCut;
    if (a < 0)
        a = 0;
    b = mean.val[0] + std.val[0] * autoLevelCut;
    m = 1. / (b - a);
    c = -m * a;

    for (y = 0, i = 0; y < ridges.rows; y++) {
        for (x = 0; x < ridges.cols; x++, i++) {
            if (r[i].mask) {
                p = (float)(m * r[i].mag + c);
                if (p > 1.0f)
                    p = 1.0f;
                else if (p < 0)
                    p = 0;

                r[i].mag = p;
                ridges.at<float>(y, x) = p;
            }
        }
    }
}

void fastBilateralFilter(cv::Mat &input, cv::Mat &output, int spatialWindow,
                         float graySigma, int order) {
    const int bins = 64;

    spatialWindow |= 1;
    graySigma *= bins;

    int i, j, o, k;
    int imW, imH, imin, imax, jmin, jmax, histH, histWB, histHWB;
    int hk = spatialWindow / 2;

    double minval, maxval, valscale;
    double sigma_2, c;
    float grayk[bins];
    float conv, norm, term;
    int line_hist[bins];
    int diff;
    int *im, *im_imin_jmin, *im_imin_jmax, *im_imax_jmin, *im_imax_jmax,
        *im_i_jmax, *im_i_j, *im_ii_j, *im_i_jj, *im_ii_jj, *im_imin_j,
        *im_imax_j, *im_i_jmin, *im_i, *im_ci, *im_ci_cj;
    int *ihist, *ihist_h, *ihist_h_1, *ihist_wh, *ihist_h_j, *ihist_h_1_j,
        *ihist_wh_j, *ihist_h_wj, *ihist_wh_wj, *ihist_end;
    int hist;
    uchar *row, *col;

    // Get value range & statistics
    cv::minMaxLoc(input, &minval, &maxval);
    valscale = (bins - 1) / (maxval - minval);

    // Init Gaussian gray kernel
    c = 1.0 / (sqrt(2.0 * PI) * graySigma);
    sigma_2 = 2.0 * graySigma * graySigma;
    for (i = 0; i < bins; i++) {
        grayk[i] = (float)(c * exp(-(double)(i * i) / sigma_2));
    }

    // Init integral image & histogram
    imW = input.cols + spatialWindow + 1;
    imH = input.rows + spatialWindow + 1;
    histH = spatialWindow + 1;
    histWB = imW * bins;
    histHWB = histH * histWB;
    jmin = imin = hk + 1;
    jmax = imW - hk - 2;
    imax = imH - hk - 2;
    im = new int[imH * imW]();
    ihist = new int[histH * histWB]();
    im_imin_jmin = im + imin * imW + jmin;
    im_imin_jmax = im + imin * imW + jmax;
    im_imax_jmin = im + imax * imW + jmin;
    im_imax_jmax = im + imax * imW + jmax;
    ihist_end = ihist + histHWB;

    // Quantize image levels to integers
    for (row = input.ptr<uchar>(), im_i_jmin = im_imin_jmin;
         im_i_jmin <= im_imax_jmin; im_i_jmin += imW, row += input.step) {
        im_i_jmax = im_i_jmin + input.cols - 1;
        for (col = row, im_i_j = im_i_jmin; im_i_j <= im_i_jmax;
             im_i_j++, col += sizeof(float)) {
            *im_i_j = (int)(((double)(*(float *)col) - minval) * valscale);
        }
    }

    // Order loop (box, triangle, approx. Gaussian) filters
    for (o = 0; o <= order; o++) {
        // Prepare image borders
        for (j = jmin, im_imin_j = im_imin_jmin, im_imax_j = im_imax_jmin;
             j <= jmax; j++, im_imin_j++, im_imax_j++) {
            for (i = 1, im_i_j = im_imin_j - imW, im_ii_j = im_imax_j + imW;
                 i < imin; i++, im_i_j -= imW, im_ii_j += imW) {
                *im_i_j = *im_imin_j;
                *im_ii_j = *im_imax_j;
            }
        }
        for (i = imin, im_i_jmin = im_imin_jmin, im_i_jmax = im_imin_jmax;
             i <= imax; i++, im_i_jmin += imW, im_i_jmax += imW) {
            for (j = 1, im_i_j = im_i_jmin - 1, im_i_jj = im_i_jmax + 1;
                 j < jmin; j++, im_i_j--, im_i_jj++) {
                *im_i_j = *im_i_jmin;
                *im_i_jj = *im_i_jmax;
            }
        }
        for (i = 1, im_i_j = im_imin_jmin - imW, im_ii_j = im_imax_jmin + imW,
            im_i_jj = im_imin_jmax - imW, im_ii_jj = im_imax_jmax + imW;
             i < imin; i++, im_i_j -= imW, im_ii_j += imW, im_i_jj -= imW,
            im_ii_jj += imW) {
            for (j = 1; j < jmin; j++) {
                *(im_i_j - j) = *im_imin_jmin;
                *(im_ii_j - j) = *im_imax_jmin;
                *(im_i_jj + j) = *im_imin_jmax;
                *(im_ii_jj + j) = *im_imax_jmax;
            }
        }

        // Init pointers
        im_i = im + imW;
        im_ci = im + (1 - hk) * imW;
        ihist_h = ihist + histWB;
        ihist_h_1 = ihist;
        ihist_wh = ihist + 2 * histWB;
        row = output.ptr<uchar>() + (1 - spatialWindow) * input.step;

        // Filter loop
        for (i = 1; i < imH - 1; i++) {
            // Clear temp histogram item
            fill(&line_hist[0], &line_hist[0] + bins, 0);

            for (j = 1, im_i_j = im_i + 1, im_ci_cj = im_ci + 1 - hk,
                ihist_h_j = ihist_h + bins, ihist_h_1_j = ihist_h_1 + bins,
                ihist_wh_j = ihist_wh + bins,
                ihist_h_wj = ihist_h + (1 - spatialWindow) * bins,
                ihist_wh_wj = ihist_wh + (1 - spatialWindow) * bins,
                col = row + (1 - spatialWindow) * sizeof(float);
                 j < imW - 1; j++, im_i_j++, im_ci_cj++, ihist_h_j += bins,
                ihist_h_1_j += bins, ihist_wh_j += bins, ihist_h_wj += bins,
                ihist_wh_wj += bins, col += sizeof(float)) {
                // Calculate integral histogram
                line_hist[*im_i_j]++;
                for (k = 0; k < bins; k++)
                    ihist_h_j[k] = line_hist[k] + ihist_h_1_j[k];

                // Calculate filtered image & hist
                if (i >= spatialWindow && j >= spatialWindow) {
                    norm = 0;
                    conv = 0;
                    for (k = 0; k < bins; k++) {
                        hist = ihist_h_j[k] - ihist_wh_j[k] - ihist_h_wj[k] +
                               ihist_wh_wj[k];
                        diff = (*im_ci_cj > k) ? *im_ci_cj - k : k - *im_ci_cj;
                        term = hist * grayk[diff];
                        conv += term * k;
                        norm += term;
                    }
                    *(float *)col = conv / norm;
                }
            }

            // Update pointers
            im_i += imW;
            im_ci += imW;
            ihist_h += histWB;
            if (ihist_h >= ihist_end)
                ihist_h -= histHWB;
            ihist_h_1 += histWB;
            if (ihist_h_1 >= ihist_end)
                ihist_h_1 -= histHWB;
            ihist_wh += histWB;
            if (ihist_wh >= ihist_end)
                ihist_wh -= histHWB;
            row += output.step;
        }
    }

    // Free mem.
    delete[] im;
    delete[] ihist;
}

void saveRidgeStructure(int s, const string tag, vector<ridge> &r) {
    if (!r.empty()) {
        string outputfn;
        outputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
        outputfn = (boost::format(outputfn) % s).str();
        int numRecords = r.size();
        cout << "Saving: " << outputfn << endl;
        ofstream f(outputfn, ios::binary);
        f.write(reinterpret_cast<char *>(&numRecords), sizeof(int));
        f.write(reinterpret_cast<char *>(r.data()), sizeof(ridge) * numRecords);
        f.close();
    }
}

vector<ridge> loadRidgeStructure(int s, const string &tag) {
    string inputfn;
    inputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
    inputfn = (boost::format(inputfn) % s).str();
    cout << "Loading: " << inputfn << endl;
    ifstream f(inputfn, ios::binary);
    int numRecords;
    f.read(reinterpret_cast<char *>(&numRecords), sizeof(int));
    vector<ridge> r(numRecords);
    if (!f.read(reinterpret_cast<char *>(r.data()), sizeof(ridge) * numRecords))
        cerr << "Error loading data!" << endl;
    f.close();
    return r;
}
