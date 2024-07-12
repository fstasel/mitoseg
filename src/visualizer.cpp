/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "visualizer.h"

void visualizeSnake(cv::Mat &visual, const double snake[][2],
                    const double scale, cv::Scalar cl, int markerSize,
                    cv::Scalar cl2) {
    int i, j, x, y, xn, yn;
    for (i = 0; i < SN_N; i++) {
        j = (i + 1) % SN_N;
        x = cvRound(scale * snake[i][0]);
        y = cvRound(scale * snake[i][1]);
        xn = cvRound(scale * snake[j][0]);
        yn = cvRound(scale * snake[j][1]);
        cv::line(visual, cv::Point(x, y), cv::Point(xn, yn), cl2);
        cv::circle(visual, cv::Point(x, y), markerSize, cl, cv::FILLED);
    }
}

Visualizer::Visualizer() {
    //
    buf.resize(MAXBUFFER);
    bufIndex.resize(MAXBUFFER);
    //
    startSlice = 0;
    numSlices = 1;
    //
    gx = gy = 1;
    gw = gh = w = h = width = height = 0;
    //
    for (int i = 0; i < MAXBUFFER; i++) {
        bufIndex[i] = -1;
    }
    nextFree = 0;
    //
    setFileNameTag("", "", "");
    //
    snake_scale = 1.0;
    //
    marker_z = -1;
    marker_scale = 1.0;
    //
    poly_scale = 1.0;
}

Visualizer::~Visualizer() { flushBuffer(); }

void Visualizer::setSnakeBuf(vector<snake25d> &buf, double scale) {
    snakeBuf = buf;
    snake_scale = scale;
}

void Visualizer::setMarkerBuf(vector<cv::Point> &buf, int z, double scale) {
    markerBuf = buf;
    marker_z = z;
    marker_scale = scale;
}

void Visualizer::clearMarkerBuf() {
    markerBuf.clear();
    marker_z = -1;
    marker_scale = 1.0;
}

void Visualizer::setPolyBuf(vector<poly25d> &buf, double scale) {
    polyBuf = buf;
    poly_scale = scale;
}

void Visualizer::setStartSlice(int s) {
    if (s >= SLICE_START && s <= SLICE_END)
        startSlice = s;
}

int Visualizer::getStartSlice() { return startSlice; }

void Visualizer::setNumSlices(int n) {
    if (n >= 1 && n <= MAXIMAGES)
        numSlices = n;
}

int Visualizer::getNumSlices() { return numSlices; }

cv::Size Visualizer::getFrameSize() { return cv::Size(width, height); }

cv::Size Visualizer::getGridSize() { return cv::Size(gw, gh); }

cv::Size Visualizer::getOutputSize() { return cv::Size(gx, gy); }

void Visualizer::setFileNameTag(const string &path, const string &tag1,
                                const string &tag2) {
    fname_tag1 = tag1;
    fname_tag2 = tag2;
    fname_path = path;
    flushBuffer();
    if (fname_tag1[0] == 0 && fname_tag2[0] == 0 && fname_path[0] == 0) {
        sliceSize = cv::Size(ROI_W, ROI_H);
    } else {
        cv::Mat buffer;
        for (int s = SLICE_START; s <= SLICE_END; s++) {
            buffer = getSlice(s);
            if (!buffer.empty()) {
                sliceSize = cv::Size(buffer.cols, buffer.rows);
                return;
            }
        }

        if (buffer.empty())
            sliceSize = cv::Size(50, 50);
    }
}

int Visualizer::getBuffer(int s, cv::Mat &im) {
    for (int i = 0; i < MAXBUFFER; i++) {
        if (bufIndex[i] == s) {
            im = buf[i];
            return 0;
        }
    }
    return -1;
}

void Visualizer::insertBuffer(int s, cv::Mat im) {
    int i, f = 0, ind = nextFree;
    for (i = 0; i < MAXBUFFER; i++) {
        if (bufIndex[i] == s) {
            ind = i;
            f = 1;
            break;
        }
    }

    buf[ind] = im.clone();
    bufIndex[ind] = s;

    if (!f)
        nextFree = (nextFree + 1) % MAXBUFFER;
}

int Visualizer::loadSliceImage(int s, cv::Mat &im) {
    cv::Mat temp;
    string destpath, fname;

    if (fname_tag1[0] == 0 && fname_tag2[0] == 0 && fname_path[0] == 0) {
        temp = loadSlice(s);
        if (!temp.empty()) {
            im = getImageROI(temp);
        }
        return temp.empty();
    } else {
        destpath = (fname_path[0] != 0) ? fname_path : DESTPATH;
        fname = (boost::format("%s%s%s%s") % destpath % fname_tag1 %
                 fname_tag2 % FNAME)
                    .str();
        fname = (boost::format(fname) % s).str();
        cout << "Loading: " << fname << endl;
        im = cv::imread(fname);
        if (im.empty()) {
            cerr << "Error loading image!" << endl;
            return 1;
        }
        if (im.channels() == 1) {
            cv::cvtColor(im, im, cv::COLOR_GRAY2BGR);
        }
        return 0;
    }
}

cv::Mat Visualizer::getSlice(int slice) {
    cv::Mat buffer;
    if (slice >= SLICE_START && slice <= SLICE_END) {
        loadSliceImage(slice, buffer);
    }
    return buffer;
}

cv::Mat Visualizer::getFrame(int slice, cv::Size size) {
    cv::Mat temp, buffer;
    if (!getBuffer(slice, buffer)) {
        if (buffer.cols == size.width && buffer.rows == size.height)
            return buffer.clone();
    }

    buffer = getSlice(slice);
    if (!buffer.empty()) {
        cv::resize(buffer, temp, size, 0, 0, cv::INTER_AREA);
        insertBuffer(slice, temp);
        return temp;
    }

    getNATemplate(buffer);
    cv::resize(buffer, temp, size, 0, 0, cv::INTER_NEAREST);
    return temp;
}

void Visualizer::getNATemplate(cv::Mat &buffer) {
    cv::Size text_size;

    buffer = cv::Mat(50, 50, CV_8UC3, cv::Scalar(0, 0, 0));
    cv::putText(buffer, "N/A", cv::Point(25, 25), cv::FONT_HERSHEY_SIMPLEX, 0.3,
                cv::Scalar(0, 255, 255));
}

void Visualizer::drawSnake(cv::Mat &buffer, snake25d &snake, int snakeSlice) {
    cv::Scalar cl;
    int marker;
    switch (snake.status) {
    case SNAKE_GROWING:
        cl = cv::Scalar(0, 255, 255);
        break;
    case SNAKE_CONVERGED:
        if (snake.isValid[snakeSlice] == SNAKE_VALID)
            cl = cv::Scalar(0, 255, 0);
        else
            cl = cv::Scalar(0, 0, 255);
        break;
    case SNAKE_MAXAREA:
    case SNAKE_MINAREA:
    case SNAKE_MAXITER:
        cl = cv::Scalar(0, 0, 255);
        break;
    default:
        cl = cv::Scalar(255, 255, 255);
    }
    marker = MAX(cvRound(5.0 * gw / width), 1);
    visualizeSnake(buffer, snake.node[snakeSlice], snake_scale * gw / width, cl,
                   marker);
}

void Visualizer::drawMarker(cv::Mat &buffer) {
    uint i;
    int x, y;
    int s = MAX(cvRound(10.0 * gw / width), 1);
    int t = MAX(s / 2, 1);
    for (i = 0; i < markerBuf.size(); i++) {
        x = cvRound(markerBuf[i].x * marker_scale * gw / width);
        y = cvRound(markerBuf[i].y * marker_scale * gw / width);
        cv::line(buffer, cv::Point(x - s, y), cv::Point(x + s, y),
                 cv::Scalar(0, 255, 255), t);
        cv::line(buffer, cv::Point(x, y - s), cv::Point(x, y + s),
                 cv::Scalar(0, 255, 255), t);
    }
}

void Visualizer::drawPoly(cv::Mat &buffer, poly25d &poly, int polySlice) {
    int t = MAX(cvRound(5.0 * gw / width), 1);
    int x1, y1, x2, y2;
    for (polygon &seq : poly.slice[polySlice]) {
        for (uint i = 0; i < seq.size(); i++) {
            x1 = cvRound(seq[i].x * poly_scale * gw / width);
            y1 = cvRound(seq[i].y * poly_scale * gw / width);
            x2 = cvRound(seq[(i + 1) % seq.size()].x * poly_scale * gw / width);
            y2 = cvRound(seq[(i + 1) % seq.size()].y * poly_scale * gw / width);
            cv::line(buffer, cv::Point(x1, y1), cv::Point(x2, y2),
                     cv::Scalar(255, 0, 0), t);
        }
    }
}

void Visualizer::drawFrame(int slice, cv::Mat &buffer) {
    uint n;
    // Draw snakes
    if (!snakeBuf.empty()) {
        int ss;
        for (n = 0; n < snakeBuf.size(); n++) {
            if (slice >= snakeBuf[n].start_z && slice <= snakeBuf[n].end_z) {
                ss = slice - snakeBuf[n].start_z;
                drawSnake(buffer, snakeBuf[n], ss);
            }
        }
    }
    // Draw markers
    if (!markerBuf.empty() && slice == marker_z) {
        drawMarker(buffer);
    }
    // Draw poly
    if (!polyBuf.empty()) {
        int ss;
        for (n = 0; n < polyBuf.size(); n++) {
            if (slice >= polyBuf[n].start_z && slice <= polyBuf[n].end_z) {
                ss = slice - polyBuf[n].start_z;
                drawPoly(buffer, polyBuf[n], ss);
            }
        }
    }
}

void Visualizer::flushBuffer() {
    int i;
    for (i = 0; i < MAXBUFFER; i++) {
        bufIndex[i] = -1;
    }
    nextFree = 0;
}

void Visualizer::update() {
    cv::Mat buffer;
    cv::Size text_size;
    int baseline;

    width = sliceSize.width;
    height = sliceSize.height;

    gx = gy = 1;
    double r = (double)WINDOW_WIDTH / WINDOW_HEIGHT;
    double ri = (double)width / height;
    while (gx * gy < numSlices) {
        double rgx = ri * (gx + 1) / gy;
        double rgy = ri * gx / (gy + 1);
        if (fabs(rgx - r) < fabs(rgy - r))
            gx++;
        else
            gy++;
    }

    w = gx * width;
    h = gy * height;
    gw = width;
    gh = height;
    if (w > WINDOW_WIDTH || h > WINDOW_HEIGHT) {
        int _h = h * WINDOW_WIDTH / w;
        if (_h <= WINDOW_HEIGHT) {
            gw = gw * WINDOW_WIDTH / w;
            gh = gh * WINDOW_WIDTH / w;
        } else {
            gw = gw * WINDOW_HEIGHT / h;
            gh = gh * WINDOW_HEIGHT / h;
        }
    }
    w = gx * gw;
    h = gy * gh;

    outputImage = cv::Mat::zeros(h, w, CV_8UC3);
    string fn;

    int x, y, px, py, i, s;
    for (i = 0, y = 0, py = 0, s = startSlice;
         y < gy && i < numSlices && s <= SLICE_END; y++, py += gh) {
        for (x = 0, px = 0; x < gx && i < numSlices && s <= SLICE_END;
             x++, px += gw, i++, s++) {
            // Frame
            buffer = getFrame(s, cv::Size(gw, gh));
            drawFrame(s, buffer);
            cv::Rect roi(px, py, gw, gh);
            buffer.copyTo(outputImage(roi));

            fn = (boost::format("%d") % (i + startSlice)).str();
            text_size = cv::getTextSize(fn, cv::FONT_HERSHEY_SIMPLEX, 0.5, 1,
                                        &baseline);
            cv::putText(outputImage, fn,
                        cv::Point(px + 2, py + text_size.height + 2),
                        cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 255, 255));
        }
    }
    // Post drawings
    for (y = 1, py = gh; y < gy; y++, py += gh) {
        cv::line(outputImage, cv::Point(0, py), cv::Point(w - 1, py),
                 cv::Scalar(0, 0, 255));
    }
    for (x = 1, px = gw; x < gx; x++, px += gw) {
        cv::line(outputImage, cv::Point(px, 0), cv::Point(px, h - 1),
                 cv::Scalar(0, 0, 255));
    }
}

void Visualizer::getSliceCoor(cv::Point p, int *s, cv::Point &c) {
    *s = -1;
    c = cv::Point(0, 0);
    if (!outputImage.empty()) {
        int x, y, t, l;
        x = p.x / gw;
        y = p.y / gh;
        l = (p.x % gw) * width / gw;
        t = (p.y % gh) * height / gh;
        *s = startSlice + y * gx + x;
        if (*s <= SLICE_END && *s < startSlice + numSlices) {
            c = cv::Point(l, t);
        } else {
            *s = -1;
            c = cv::Point(0, 0);
        }
    }
}
