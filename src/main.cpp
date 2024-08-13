/*
 *  SPDX-FileCopyrightText: 2024 Faris Serdar Taşel <fst@cankaya.edu.tr>
 *  SPDX-FileCopyrightText: 2024 Efe Çiftci <efeciftci@cankaya.edu.tr>
 *
 *  SPDX-License-Identifier: GPL-3.0-or-later
 */

#define VERSION "1.0"

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <opencv2/opencv.hpp>

#include "phase1_main.h"
#include "phase2_main.h"
#include "phase3_main.h"
#include "settings_parser.h"
#include "sys/resource.h"
using namespace std;
using namespace mitoseg_settings;
namespace po = boost::program_options;

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

int main(int argc, char **argv) {
    po::options_description options("Usage");
    options.add_options()
        ("help", "display this help and exit")
        ("zrange", po::value<vector<int>>()->multitoken(), "z-range of slices")
        ("psize", po::value<double>(), "pixel size as nm/px")
        ("roi", po::value<vector<int>>()->multitoken(), "region of interest")
        ("src", po::value<string>(), "source directory")
        ("dst", po::value<string>(), "destination directory")
        ("phase", po::value<int>(), "apply a specific phase only")
        ("valid", po::value<double>(), "validity threshold between 0 - 1")
        ("thick", po::value<string>(), "snake z-thickness")
        ("cores", po::value<int>(), "# of cpu cores to utilize")
        ("settingsFile", po::value<string>(), "path of settings file")
        ("pattern", po::value<string>(), "filename pattern");

    po::positional_options_description p;
    p.add("pattern", -1);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(options)
                      .positional(p)
                      .run(),
                  vm);
    } catch (po::error &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    po::notify(vm);

    if (vm.count("help") || argc < 2) {
        cout
            << "MitoSeg v" << VERSION << endl
            << "  Faris Serdar Tasel <fst@cankaya.edu.tr" << endl
            << "  Efe Ciftci <efeciftci@cankaya.edu.tr" << endl
            << endl
            << "Based on the paper:" << endl
            << "Tasel, S.F., Mumcuoglu, E.U., Hassanpour, R.Z., Perkins, G., "
               "2016."
            << endl
            << "A validated active contour method driven by parabolic arc model"
            << endl
            << "for detection and segmentation of mitochondria." << endl
            << "J. Struct. Biol. 194, 253–271. doi:10.1016/j.jsb.2016.03.002"
            << endl
            << endl
            << "Options:" << endl
            << "  --help\t\tdisplay this help and exit" << endl
            << "  --zrange START END\tz-range of slices" << endl
            << "  --psize SIZE\t\tpixel size as nm/px" << endl
            << "  --roi X Y W H\t\tregion of interest" << endl
            << "  --src PATH\t\tsource directory" << endl
            << "  --dst PATH\t\tdestination directory" << endl
            << "  --phase PHASE\t\tapply a specific phase only (1, 2, or 3)"
            << endl
            << "  --valid VALIDITY\tvalidity threshold between 0 - 1 (default: "
               "0.75)"
            << endl
            << "  --thick THICKNESS\tsnake z-thickness (default: 20)" << endl
            << "  --cores CORES\t\t# of cpu cores to utilize" << endl
            << "  --settingsFile PATH\tpath of custom settings file" << endl
            << endl
            << "The options --zrange, --psize and <filename pattern> are "
               "mandatory."
            << endl
            << endl
            << "Example:" << endl
            << "  mitoseg --zrange 30 100 --psize 2.0 dataset_slice%04d.tif"
            << endl
            << "    Process files dataset_slice0030.tif ... "
               "dataset_slice0100.tif"
            << endl
            << "    assuming that pixel size is 2.0nm" << endl
            << endl
            << "  mitoseg --zrange 40 120 --psize 1.1 mito%d.bmp" << endl
            << "    Process files mito40.bmp ... mito120.bmp" << endl
            << "    assuming that pixel size is 1.1nm" << endl
            << endl
            << "E-mail to fst@cankaya.edu.tr for more information." << endl;
        return 0;
    }

    if (vm.count("settingsFile"))
        SETTINGS_PATH = vm["settingsFile"].as<string>();
    else
        SETTINGS_PATH = "settings.yaml";
    settingsSet = loadSettings();

    if (vm.count("zrange")) {
        SLICE_START = vm["zrange"].as<vector<int>>()[0];
        SLICE_END = vm["zrange"].as<vector<int>>()[1];
        zrangeSet = true;
        if (SLICE_START < 0 || SLICE_END < 0 || SLICE_START > SLICE_END) {
            cerr << "Error: Invalid parameter for --zrange" << endl;
            return 1;
        }
    }

    if (vm.count("psize")) {
        RESOLUTION = vm["psize"].as<double>();
        psizeSet = true;
        if (RESOLUTION <= 0) {
            cerr << "Error: Invalid parameter for --psize" << endl;
            return 1;
        }
    }

    if (vm.count("roi")) {
        ROI_X = vm["roi"].as<vector<int>>()[0];
        ROI_Y = vm["roi"].as<vector<int>>()[1];
        ROI_W = vm["roi"].as<vector<int>>()[2];
        ROI_H = vm["roi"].as<vector<int>>()[3];
        roiSet = true;
        if (ROI_X < 0 || ROI_Y < 0 || ROI_W <= 0 || ROI_H <= 0) {
            cerr << "Error: Invalid parameter for --roi" << endl;
            return 1;
        }
    }

    if (vm.count("src")) {
        SOURCEPATH = vm["src"].as<string>();
        srcSet = true;
    }

    if (vm.count("dst")) {
        DESTPATH = vm["dst"].as<string>();
        dstSet = true;
    }

    if (vm.count("phase")) {
        phaseNum = vm["phase"].as<int>();
        phaseSet = true;
        if (phaseNum < 1 || phaseNum > 3) {
            cerr << "Error: Phase # must be between 1 and 3." << endl;
            return 1;
        }
    }

    if (vm.count("valid")) {
        POLY_VALIDITY = vm["valid"].as<double>();
        validSet = true;
        if (POLY_VALIDITY < 0 || POLY_VALIDITY > 1) {
            cerr << "Error: Invalid parameter for --valid (must be between 0 - "
                    "1)"
                 << endl;
            return 1;
        }
    }

    if (vm.count("thick")) {
        if (vm["thick"].as<string>() == "full") {
            thickSet = true;
        } else {
            sn25d_t = stoi(vm["thick"].as<string>());
            if (sn25d_t < 5 || sn25d_t > 500) {
                cerr << "Error: Invalid parameter for --thick (must be between "
                        "5 - 500)"
                     << endl;
                return 1;
            }
        }
    }

    if (vm.count("cores")) {
        numCores = vm["cores"].as<int>();
        if (numCores < 1 || numCores > 256) {
            cerr << "Error: # of cores must be between 1 and 256." << endl;
            return 1;
        }
    }

    if (vm.count("pattern")) {
        FNAME = vm["pattern"].as<string>();
        fnameSet = true;
    }

    if (!(zrangeSet && psizeSet && fnameSet)) {
        cout << "The options --zrange, --psize and <filename pattern> are "
                "mandatory input."
             << endl;
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
        return 1;
    }

    if (settingsSet)
        cout << "Settings: loaded from " << SETTINGS_PATH << endl;
    else
        cout << "Settings: using defaults" << endl;

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
