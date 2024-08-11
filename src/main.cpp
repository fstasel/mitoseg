/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#define VERSION "1.0"

#include <boost/format.hpp>
#include <opencv2/opencv.hpp>

#include "phase1_main.h"
#include "phase2_main.h"
#include "phase3_main.h"
#include "settings_parser.h"
#include "sys/resource.h"
using namespace std;
using namespace mitoseg_settings;

bool zrangeSet = false;
bool psizeSet = false;
bool roiSet = false;
bool srcSet = false;
bool dstSet = false;
bool phaseSet = false;
bool fnameSet = false;
bool validSet = false;
bool thickSet = false;
bool settingsSet = false;
int phaseNum = 1;
int thickParam = 0;

void displayHelp() {
    cout << "MitoSeg v" << VERSION << endl;
    cout << "\tFaris Serdar Tasel <fst@cankaya.edu.tr>" << endl;
    cout << "\tEfe Ciftci <efeciftci@cankaya.edu.tr>" << endl << endl;
    cout << "Based on the paper:" << endl;
    cout << "Tasel, S.F., Mumcuoglu, E.U., Hassanpour, R.Z., Perkins, G., 2016."
         << endl
         << "A validated active contour method driven by parabolic arc model"
         << endl
         << "for detection and segmentation of mitochondria." << endl
         << "J. Struct. Biol. 194, 253–271. doi:10.1016/j.jsb.2016.03.002"
         << endl
         << endl;
    cout << "Usage:" << endl;
    cout << "mitoseg -options <parameters> <filename pattern>" << endl << endl;
    cout << "Options:" << endl;
    cout << "\t-zrange <start slice #> <end slice #>\tSpecify z-range" << endl;
    cout << "\t-psize <pixel size>\t\t\tSpecify pixel size as nm/px" << endl;
    cout << "\t-roi <left> <top> <width> <height>\tSpecify region of interest"
         << endl;
    cout << "\t-src <source directory>\t\t\tSpecify source directory" << endl;
    cout << "\t-dst <destination directory>\t\tSpecify destination directory"
         << endl;
    cout << "\t-phase <phase #>\t\t\tApply phase <phase #> only." << endl
         << "\t\t\t\t\t\t(must be 1, 2 or 3)" << endl;
    cout << "\t-valid <validity threshold>\t\tSpecify validity threshold"
         << endl
         << "\t\t\t\t\t\tbetween 0-1 (default: 0.75)" << endl;
    cout << "\t-thick <z-thickness>\t\t\tSpecify snake z-thickness" << endl
         << "\t\t\t\t\t\tbetween 5-500 (default: 20)" << endl
         << "\t\t\t\t\t\tset 'full' to use z-range" << endl;
    cout << "\t-cores <# of cpu cores>\t\t\tSet # of cpu cores to utilize"
         << endl
         << "\t\t\t\t\t\t(default: 1)" << endl;
    cout << "\t-settingsFile <path>\t\t\tPath of settings file" << endl
         << endl
         << endl;
    cout << "The options -zrange, -psize and <filename pattern> are mandatory "
            "input."
         << endl
         << endl;
    cout << "Example:" << endl;
    cout << "\tmitoseg -zrange 30 100 -psize 2.0 dataset_slice%04d.tif" << endl;
    cout << "\t\tProcess files dataset_slice0030.tif ... dataset_slice0100.tif"
         << endl;
    cout << "\t\tassuming that pixel size is 2.0nm" << endl;
    cout << "\tmitoseg -zrange 40 120 -psize 1.1 mito%d.bmp" << endl;
    cout << "\t\tProcess files mito40.bmp ... mito120.bmp" << endl;
    cout << "\t\tassuming that pixel size is 1.1nm" << endl << endl;
    cout << "E-mail to fst@cankaya.edu.tr for more information." << endl
         << endl;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        displayHelp();
        return 0;
    }

    int c = 1;
    while (c < argc) {
        if (strcmp(argv[c], "-settingsFile") == 0) {
            if (c + 1 < argc) {
                SETTINGS_PATH = argv[++c];
                settingsSet = true;
            } else {
                cerr << "Error: Insufficient parameter for "
                        "-settingsFile"
                     << endl
                     << "Use mitoseg with no parameters for help" << endl
                     << endl;
                return 1;
            }
        }
        c++;
    }

    if (!settingsSet)
        SETTINGS_PATH = "settings.yaml";
    if (loadSettings())
        cout << "Settings: loaded from " << SETTINGS_PATH << endl;
    else
        cout << "Settings: using defaults" << endl;

    c = 1;
    while (c < argc) {
        if (strcmp(argv[c], "-zrange") == 0) {
            if (c + 2 < argc) {
                SLICE_START = atoi(argv[++c]);
                SLICE_END = atoi(argv[++c]);
                zrangeSet = true;
                if (SLICE_START < 0 || SLICE_END < 0 ||
                    SLICE_START > SLICE_END) {
                    cerr << "Error: Invalid parameter for "
                            "-zrange"
                         << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -zrange" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-psize") == 0) {
            if (c + 1 < argc) {
                RESOLUTION = atof(argv[++c]);
                psizeSet = true;
                if (RESOLUTION <= 0) {
                    cerr << "Error: Invalid parameter for "
                            "-psize"
                         << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -psize" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-roi") == 0) {
            if (c + 4 < argc) {
                ROI_X = atoi(argv[++c]);
                ROI_Y = atoi(argv[++c]);
                ROI_W = atoi(argv[++c]);
                ROI_H = atoi(argv[++c]);
                roiSet = true;
                if (ROI_X < 0 || ROI_Y < 0 || ROI_W <= 0 || ROI_H <= 0) {
                    cerr << "Error: Invalid parameter for -roi" << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -roi" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-src") == 0) {
            if (c + 1 < argc) {
                SOURCEPATH = argv[++c];
                srcSet = true;
            } else {
                cerr << "Error: Insufficient parameter for -src" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-dst") == 0) {
            if (c + 1 < argc) {
                DESTPATH = argv[++c];
                dstSet = true;
            } else {
                cerr << "Error: Insufficient parameter for -dst" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-phase") == 0) {
            if (c + 1 < argc) {
                phaseNum = atoi(argv[++c]);
                phaseSet = true;
                if (phaseNum < 1 || phaseNum > 3) {
                    cerr << "Error: Phase # must be between 1 "
                            "and 3."
                         << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -phase" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-valid") == 0) {
            if (c + 1 < argc) {
                POLY_VALIDITY = atof(argv[++c]);
                validSet = true;
                if (POLY_VALIDITY < 0 || POLY_VALIDITY > 1) {
                    cerr << "Error: Invalid parameter for "
                            "-valid (must be "
                            "between 0-1)"
                         << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -valid" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-thick") == 0) {
            if (c + 1 < argc) {
                if (strcmp(argv[c + 1], "full") == 0) {
                    thickSet = true;
                } else {
                    sn25d_t = atoi(argv[c + 1]);
                    if (sn25d_t < 5 || sn25d_t > 500) {
                        cerr << "Error: Invalid parameter "
                                "for -thick (must be "
                                "between 5-500)"
                             << endl;
                        cerr << "Use mitoseg with no "
                                "parameter for help"
                             << endl
                             << endl;
                        return 1;
                    }
                }
                c++;
            } else {
                cerr << "Error: Insufficient parameter for -thick" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-cores") == 0) {
            if (c + 1 < argc) {
                numCores = atoi(argv[++c]);
                if (numCores < 1 || numCores > 256) {
                    cerr << "Error: # of cores must be "
                            "between 1 and 256."
                         << endl;
                    cerr << "Use mitoseg with no parameter "
                            "for help"
                         << endl
                         << endl;
                    return 1;
                }
            } else {
                cerr << "Error: Insufficient parameter for -cores" << endl;
                cerr << "Use mitoseg with no parameter for help" << endl
                     << endl;
                return 1;
            }
        } else if (strcmp(argv[c], "-settingsFile") == 0) {
            // settingsFile is already handled above
            c++;
        } else if (argv[c][0] == '-') {
            cerr << "Error: Unknown option " << argv[c] << endl;
            cerr << "Use mitoseg with no parameter for help" << endl << endl;
            return 1;
        } else {
            FNAME = argv[c];
            fnameSet = true;
        }
        c++;
    }

    if (!(zrangeSet && psizeSet && fnameSet)) {
        cout << "The options -zrange, -psize and <filename pattern> are "
                "mandatory input."
             << endl;
        cout << "Use mitoseg with no parameter for help" << endl << endl;
        return 1;
    }
    if (thickSet) {
        sn25d_t = SLICE_END - SLICE_START + 1;
        if (sn25d_t < 5)
            sn25d_t = 5;
    }
    if (SLICE_END - SLICE_START + 1 < sn25d_t) {
        cerr << "Error: Z-range must contain at least " << sn25d_t << " slices."
             << endl;
        cerr << "Use mitoseg with no parameter for help" << endl << endl;
        return 1;
    }
    cout << boost::format("Using pixel size: %.2f nm") % RESOLUTION << endl;
    if (!roiSet) {
        cout << "ROI was not set. Enabling Auto-ROI..." << endl;
        ROI_X = ROI_Y = ROI_W = ROI_H = -1;
    } else
        cout << boost::format("Using ROI: %d %d %d %d") % ROI_X % ROI_Y %
                    ROI_W % ROI_H
             << endl;
    cout << "Using snake z-thickness = " << sn25d_t << endl;
    cout << boost::format("Using validity threshold = %.2f") % POLY_VALIDITY
         << endl;
    if (!srcSet)
        SOURCEPATH = "./";
    else if (SOURCEPATH.back() != '/')
        SOURCEPATH += "/";
    cout << "Source path: " << SOURCEPATH << endl;
    if (!dstSet)
        DESTPATH = "./";
    else if (DESTPATH.back() != '/')
        DESTPATH += "/";
    cout << "Destination path: " << DESTPATH << endl;
    cout << "Utilized cpu cores: " << numCores << endl;
    cout << "Executing: ";
    if (phaseSet)
        cout << "Phase #" << phaseNum << endl;
    else
        cout << "All phases" << endl;
    cout << "Processing from ";
    cout << boost::format(FNAME) % SLICE_START;
    cout << boost::format(" (slice #:%d) to ") % SLICE_START;
    cout << boost::format(FNAME) % SLICE_END;
    cout << boost::format(" (slice #:%d)") % SLICE_END << endl;

    if (phaseSet) {
        switch (phaseNum) {
        case 1:
            startPhase1();
            break;
        case 2:
            startPhase2();
            break;
        case 3:
            startPhase3();
            break;
        default:
            break;
        }
    } else {
        startPhase1();
        startPhase2();
        startPhase3();
    }
    cout << ">>>> Completed!" << endl << endl;

    return 0;
}
