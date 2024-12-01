/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "datasets.h"

string FNAME;
string DESTPATH;
string SOURCEPATH;
int SLICE_START;
int SLICE_END;
int ROI_X;
int ROI_Y;
int ROI_W;
int ROI_H;
double RESOLUTION;

cv::Mat loadSlice(int s) {
    string inputfn;
    inputfn = (boost::format("%s%s") % SOURCEPATH % FNAME).str();
    inputfn = (boost::format(inputfn) % s).str();

    cout << "Loading: " << inputfn << endl;
    cv::Mat im = cv::imread(inputfn, cv::IMREAD_COLOR);

    if (im.empty()) {
        cerr << "Error loading image!" << endl;
        return cv::Mat();
    }
    return im;
}

void saveImage(int s, const string &tag, const string &tag2,
               const cv::Mat &im) {
    string outputfn;
    outputfn =
        (boost::format("%s%s%s%s") % DESTPATH % tag % tag2 % FNAME).str();
    outputfn = (boost::format(outputfn) % s).str();
    if (!im.empty()) {
        cout << "Saving: " << outputfn << endl;
        cv::imwrite(outputfn, im);
    }
}

void saveImageDat(int s, const string &tag, const cv::Mat &im) {
    string outputfn;
    outputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
    outputfn = (boost::format(outputfn) % s).str();
    if (!im.empty()) {
        cout << "Saving: " << outputfn << endl;
        ofstream outfile(outputfn, ios::binary);
        if (outfile.is_open()) {
            outfile.write(reinterpret_cast<const char *>(&im.dims),
                          sizeof(int));
            outfile.write(reinterpret_cast<const char *>(&im.rows),
                          sizeof(int));
            outfile.write(reinterpret_cast<const char *>(&im.cols),
                          sizeof(int));
            outfile.write(reinterpret_cast<const char *>(im.data),
                          im.total() * im.elemSize());
            outfile.close();
        }
    }
}

cv::Mat loadImageDat(int s, const string &tag) {
    int dims, sizes[2];
    string inputfn;
    inputfn = (boost::format("%s%s%s.dat") % DESTPATH % tag % FNAME).str();
    inputfn = (boost::format(inputfn) % s).str();

    cout << "Loading: " << inputfn << endl;
    ifstream infile(inputfn, ios::binary);
    infile.read(reinterpret_cast<char *>(&dims), sizeof(int));
    infile.read(reinterpret_cast<char *>(&sizes[0]), sizeof(int));
    infile.read(reinterpret_cast<char *>(&sizes[1]), sizeof(int));
    cv::Mat im(dims, sizes, CV_32F);
    infile.read(reinterpret_cast<char *>(im.data), im.total() * im.elemSize());
    infile.close();
    return im;
}

double operator~(const cv::Vec3b &a) {
    double d;
    d = a.val[0] * a.val[0];
    d += a.val[1] * a.val[1];
    d += a.val[2] * a.val[2];
    return d;
}

double getDev(const cv::Mat &input) {
    double d = 0;
    for (int y = 1; y < input.rows - 1; y++) {
        for (int x = 1; x < input.cols - 1; x++) {
            d += ~(0.5 * input.at<cv::Vec3b>(y, x) -
                   0.125 * input.at<cv::Vec3b>(y - 1, x) -
                   0.125 * input.at<cv::Vec3b>(y + 1, x) -
                   0.125 * input.at<cv::Vec3b>(y, x - 1) -
                   0.125 * input.at<cv::Vec3b>(y, x + 1));
        }
    }
    return d / (input.rows * input.cols);
}

bool checkRow(const cv::Mat &input, const int y, const double &dev) {
    double d = 0;
    int yp, ym;
    yp = ym = y;
    if (y > 0)
        ym--;
    if (y < input.rows - 1)
        yp++;
    for (int x = 0; x < input.cols; x++)
        d += ~(0.5 * input.at<cv::Vec3b>(y, x) -
               0.25 * input.at<cv::Vec3b>(yp, x) -
               0.25 * input.at<cv::Vec3b>(ym, x));
    d /= input.cols;
    return 6.0 * d < dev;
}

bool checkColumn(const cv::Mat &input, const int x, const double &dev) {
    double d = 0;
    int xp, xm;
    xp = xm = x;
    if (x > 0)
        xm--;
    if (x < input.cols - 1)
        xp++;
    for (int y = 0; y < input.rows; y++)
        d += ~(0.5 * input.at<cv::Vec3b>(y, x) -
               0.25 * input.at<cv::Vec3b>(y, xp) -
               0.25 * input.at<cv::Vec3b>(y, xm));
    d /= input.rows;
    return 6.0 * d < dev;
}

void autoROI(const cv::Mat &input) {
    double dev = getDev(input);
    ROI_X = ROI_Y = 0;
    ROI_W = input.cols;
    ROI_H = input.rows;
    while (ROI_H > 1 && checkRow(input, ROI_Y, dev)) {
        ROI_Y++;
        ROI_H--;
    }
    while (ROI_H > 1 && checkRow(input, ROI_Y + ROI_H - 1, dev)) {
        ROI_H--;
    }
    while (ROI_W > 1 && checkColumn(input, ROI_X, dev)) {
        ROI_X++;
        ROI_W--;
    }
    while (ROI_W > 1 && checkColumn(input, ROI_X + ROI_W - 1, dev)) {
        ROI_W--;
    }
    ROI_X += 5;
    ROI_Y += 5;
    ROI_W -= 10;
    ROI_H -= 10;
}

cv::Mat getImageROI(cv::Mat &input) {
    if (ROI_X < 0) {
        autoROI(input);
        cout << boost::format("Using ROI: %d %d %d %d") % ROI_X % ROI_Y %
                    ROI_W % ROI_H
             << endl;
    }
    cv::Rect roiRect(ROI_X, ROI_Y, ROI_W, ROI_H);
    return input(roiRect).clone();
}
