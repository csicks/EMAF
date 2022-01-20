/***************************************************************************
 *
 * Authors:    Yuxuan Chen
 *
 * Copyright (C) 2021 Pattern Recognition and Bioinformatics Group, Shanghai Jiao Tong University
 *
 * Licensed under the GNU General Public License v3.0 (see LICENSE for details)
 *
 * All comments concerning this program package may be sent to e-mail address 'yxchen11@sjtu.edu.cn'
 ***************************************************************************/

#include <iostream>
#include <map>
#include "mrc.h"
#include "image2D.h"
#include "align.h"
#include <ctime>
#include "matrix.h"
#include "xmippXMD.h"
#include "noise.h"
#include "alignEvaluation.h"

using namespace std;

#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

bool getParams(int argc, char *argv[], map<string, string> &params)
{
    bool helpFlag = false;
    for (int i = 1; i < argc; i += 1)
    {
        string str = argv[i];
        if (str == "-h")
        {
            cout << BOLDRED << "USAGE: " << RESET << endl;
            cout << "\tAlign a set of images." << endl;

            cout << BOLDRED << "OPTIONS: " << RESET << endl;
            cout << BLUE << "[-h]" << RESET << endl;
            cout << "\t Show commandline instructions." << endl;
            cout << BLUE << "[-i <string>]" << RESET << endl;
            cout << "\t Path of images to be aligned." << endl;
            cout << BLUE << "[-o <string>]" << RESET << endl;
            cout << "\t Path of output aligned images." << endl;
            cout << BLUE << R"([-m <"LM" or "MD" or "DA">])" << RESET << endl;
            cout << "\t Method used for alignment. LM: Local Maximum. MD: Main Direction. DA: Direct Alignment" << endl;
            cout << BLUE << R"([-e <"LM" or "MD">])" << RESET << endl;
            cout
                << "\t Estimation of parameters of alignment method using all data. LM: Local Maximum. MD: Main Direction. Only mandatory parameters (wl for Main Direction method and wh for Local Maximum method) will be estimated."
                << endl;
            cout << GREEN << R"([-r <string: "XMIPP" or "RFAA" or path>])" << RESET << endl;
            cout
                << "\t Reference image for alignment. XMIPP for XMIPP's way of generating reference. RFAA for SPIDER's way of generating reference. Path for given path of reference image. If none, XMIPP's way will be used."
                << endl;
            cout << GREEN << "[-wh <integer>]" << RESET << endl;
            cout << "\t Radius for high-pass filter. Required by Local Maximum method and optional for Main Direction method. If none, automatic estimation using 1/10 data will be executed for Local Maximum method and value set to zero for Main Direction method." << endl;
            cout << GREEN << "[-wl <integer>]" << RESET << endl;
            cout << "\t Radius for low-pass filter. Required by Main Direction method. If none, automatic estimation using 1/10 data will be executed." << endl;
            cout << GREEN << "[-th <float>]" << RESET << endl;
            cout << "\t Threshold for power spectrum cutting. Used by Main Direction method. If none, default value will be used." << endl;
            cout << GREEN << R"([-f <"True" or "False">])" << RESET << endl;
            cout
                << "\t Whether to use fine search. Note fine search takes much more time so that it should only be used when normal ways doesn't work. If none, it will be set to False."
                << endl;
            cout << GREEN << R"([-x <"True" or "False">])" << RESET << endl;
            cout
                << "\t Whether to output XMIPP middle file. Set it to true if you would like to use this alignment method in NCEM, SPREAD or other program that makes use of XMIPP libraries. If none, it will be set to False."
                << endl;
            cout << GREEN << R"([-u <"True" or "False">])" << RESET << endl;
            cout
                << "\t Whether to output update reference during progress. Note it should be set to True if XMIPP's way of generating reference is used. If none, it will be set to True."
                << endl;

            cout << BOLDRED << "EXAMPLES: " << RESET << endl;
            cout << CYAN << "Estimate width: " << RESET << endl;
            cout << "\t ./EMAF -i /path/to/input -o /path/to/output -e LM" << endl;
            cout << CYAN << "Align using Local Maximum: " << RESET << endl;
            cout << "\t ./EMAF -i /path/to/input -o /path/to/output -m LM -wh 6" << endl;
            cout << CYAN << "Align using Main Direction: " << RESET << endl;
            cout << "\t ./EMAF -i /path/to/input -o /path/to/output -m MD -wl 10" << endl;
            helpFlag = true;
            break;
        }
        else if (str == "-i")
        {
            params.insert(pair<string, string>("-i", argv[i + 1]));
            i += 1;
        }
        else if (str == "-o")
        {
            params.insert(pair<string, string>("-o", argv[i + 1]));
            i += 1;
        }
        else if (str == "-m")
        {
            params.insert(pair<string, string>("-m", argv[i + 1]));
            i += 1;
        }
        else if (str == "-e")
        {
            params.insert(pair<string, string>("-e", argv[i + 1]));
            i += 1;
        }
        else if (str == "-r")
        {
            params.insert(pair<string, string>("-r", argv[i + 1]));
            i += 1;
        }
        else if (str == "-wl")
        {
            params.insert(pair<string, string>("-wl", argv[i + 1]));
            i += 1;
        }
        else if (str == "-wh")
        {
            params.insert(pair<string, string>("-wh", argv[i + 1]));
            i += 1;
        }
        else if (str == "-th")
        {
            params.insert(pair<string, string>("-th", argv[i + 1]));
            i += 1;
        }
        else if (str == "-f")
        {
            params.insert(pair<string, string>("-f", argv[i + 1]));
            i += 1;
        }
        else if (str == "-x")
        {
            params.insert(pair<string, string>("-f", argv[i + 1]));
            i += 1;
        }
        else if (str == "-u")
        {
            params.insert(pair<string, string>("-u", argv[i + 1]));
            i += 1;
        }
        else
            throw baseException("Error: Unsupported input type! Use -h to see commandline instructions.");
    }
    if (!helpFlag)
    {
        if (params.count("-i") == 0)
            throw baseException("Error: Input path is in need! Use -h to see commandline instructions.");
        if (params.count("-o") == 0)
            throw baseException("Error: Output path is in need!");
        if (params.count("-m") == 0 && params.count("-e") == 0)
            throw baseException("Error: One task should be specified, estimation or alignment! Use -h to see commandline instructions.");
        if (params.count("-m") != 0 && params.count("-e") != 0)
            throw baseException("Error: Only one task should be specified, estimation or alignment! Use -h to see commandline instructions.");
    }
    return helpFlag;
}

/// Please see ./EMAF -h for input format instructions.
/// Though some xmd files are supported, it is strongly suggested to convert input data into MRC format.
/// When estimation of filter's parameters is not accurate, manually fine tuning might help.
/// The source code is compiled into RELEASE version by default as specified in CMakeLists.txt. You may change it if DEBUG purpose is needed.
/// If you would like to look into the code and want to run the program with less RAM, you may change "stackReal" in code to "fileStack" which would make use of the disk and greatly reduce the demand for RAM. Such change also helps with the big data size problem.
/// Evaluation functions are provided in 'evaluation' folder with proper comments.
int main(int argc, char *argv[])
{
    map<string, string> params;
    /// As shown below, you may comment the below line and directly set params to avoid commandline inputs.
    /**
     * bool helpFlag = false;
     * params.insert(pair<string, string>("-i", "/path/to/input"));
     * params.insert(pair<string, string>("-o", "/path/to/output"));
     * params.insert(pair<string, string>("-m", "MD"));
     * **/
    bool helpFlag = getParams(argc, argv, params);
    if (helpFlag)
        return 0;
    else
    {
        bool XMIPP;
        if (params.count("-m") != 0)
        {
            if (params.count("-x") == 0)
                XMIPP = false;
            else if (params["-x"] == "True")
                XMIPP = true;
            else if (params["-x"] == "False")
                XMIPP = false;
            else
                throw baseException("Error: Unsupported XMIPP flag! Use -h to see commandline instructions.");
        }

        string pathIn = params["-i"];
        string pathOut = params["-o"];
        cout << "Input path: " << pathIn << endl;
        cout << "Output path: " << pathOut << endl;

        std::istringstream sin(pathIn);
        std::vector<std::string> fields;
        if (pathIn[0] == '.')
            fields.emplace_back(".");
        std::string field;
        while (getline(sin, field, '.'))
        {
            fields.emplace_back(field);
        }

        std::istringstream sin2(pathOut);
        std::vector<std::string> fields2;
        if (pathOut[0] == '.')
            fields2.emplace_back(".");
        std::string field2;
        while (getline(sin2, field2, '.'))
        {
            fields2.emplace_back(field2);
        }

        std::string name;
        for (int i = 0; i < fields2.size() - 1; ++i)
            name += fields2[i];

        stackReal<double> data;
        vector<string> info;
        /// "xmd" file for certain XMIPP format. Again, MRC file is strongly recommended.
        if (fields[fields.size() - 1] == "xmd")
        {
            data = readXMD(pathIn, info);
        }
        else if (fields[fields.size() - 1] == "mrcs" || fields[fields.size() - 1] == "mrc")
        {
            mrcFile mrc = mrcFile();
            mrc.read(pathIn);
            data = mrc.data;
            if (XMIPP)
                generateInfo(pathIn, data.shape[0], info);
        }
        else
        {
            throw baseException("Error: Unsupported file format! Only MRC or XMD files are supported.");
        }

        if (params.count("-m") == 0)
        {
            string nameMethod = params["-e"];
            clock_t st, et;
            st = clock();
            if (nameMethod == "LM")
            {
                int width = circleEstimatePeak(data);
                cout << "Estimated high-pass width for Local Maximum method: " << width << endl;
                return 0;
            }
            else if (nameMethod == "MD")
            {
                int width = circleEstimateShape(data);
                cout << "Estimated low-pass width for Main Direction method: " << width << endl;
                return 0;
            }
            else
                throw baseException("Error: Unsupported method type! Use -h to see commandline instructions.");
            et = clock();
            std::cout << "Estimation time:" << (double)(et - st) / CLOCKS_PER_SEC << std::endl;
        }
        else
        {
            string nameMethod = params["-m"];
            clock_t st, et;
            st = clock();
            imageReal<double> ref;
            bool fineFlag;
            if (params.count("-f") == 0)
                fineFlag = false;
            else if (params["-f"] == "False")
                fineFlag = false;
            else if (params["-f"] == "True")
                fineFlag = true;
            else
                throw baseException("Error: Unsupported fine search flag! Use -h to see commandline instructions.");
            bool updateFlag;
            if (params.count("-u") == 0)
                updateFlag = true;
            else if (params["-u"] == "False")
                updateFlag = false;
            else if (params["-u"] == "True")
                updateFlag = true;
            else
                throw baseException("Error: Unsupported update flag! Use -h to see commandline instructions.");
            if (nameMethod == "LM")
            {
                int width;
                if (params.count("-wh") == 0)
                {
                    stackReal<double> temp = data.pieces(0, int(data.shape[0] / 10));
                    width = circleEstimatePeak(temp, fineFlag);
                    cout << "Estimated high-pass width: " << width << endl;
                }
                else
                {
                    width = std::stoi(params["-wh"]);
                    cout << "Input high-pass width: " << width << endl;
                }

                if (params.count("-r") == 0)
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "XMIPP")
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "RFAA")
                {
                    ref = reference(data, alignPeak, width, fineFlag);
                }
                else
                {
                    string pathRef = params["-r"];
                    mrcFile mrcRef = mrcFile();
                    mrcRef.read(pathIn);
                    ref = mrcRef.data.pieceGet(0);
                }

                stackReal<double> stkBackup;
                if (XMIPP)
                    stkBackup = data;

                for (int i = 0; i < 2; ++i)
                    refinement(data, ref, alignPeak, width, fineFlag, updateFlag);

                ref = data.stackAvg();

                if (XMIPP)
                    outInfo(stkBackup, ref, alignPeak, info, name, width, fineFlag);
            }
            else if (nameMethod == "MD")
            {
                int outWidth, inWidth;
                double th;

                if (params.count("-wh") == 0)
                {
                    inWidth = -1;
                }
                else
                {
                    inWidth = std::stoi(params["-wh"]);
                    cout << "Input high-pass width: " << inWidth << endl;
                }

                if (params.count("-th") == 0)
                {
                    th = -1;
                }
                else
                {
                    th = std::stof(params["-th"]);
                    cout << "Input threshold: " << th << endl;
                }

                if (params.count("-wl") == 0)
                {
                    stackReal<double> temp = data.pieces(0, int(data.shape[0] / 10));
                    outWidth = circleEstimateShape(temp, inWidth, th, fineFlag);
                    cout << "Estimated low-pass width: " << outWidth << endl;
                }
                else
                {
                    outWidth = std::stoi(params["-wl"]);
                    cout << "Input low-pass width: " << outWidth << endl;
                }

                if (params.count("-r") == 0)
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "XMIPP")
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "RFAA")
                {
                    ref = reference(data, alignShape, outWidth, inWidth, th, fineFlag);
                }
                else
                {
                    string pathRef = params["-r"];
                    mrcFile mrcRef = mrcFile();
                    mrcRef.read(pathIn);
                    ref = mrcRef.data.pieceGet(0);
                }

                stackReal<double> stkBackup;
                if (XMIPP)
                    stkBackup = data;

                for (int i = 0; i < 2; ++i)
                    refinement(data, ref, alignShape, outWidth, inWidth, th, fineFlag, updateFlag);

                ref = data.stackAvg();

                if (XMIPP)
                    outInfo(stkBackup, ref, alignShape, info, name, outWidth, inWidth, th, fineFlag);
            }
            else if (nameMethod == "DA")
            {
                if (params.count("-r") == 0)
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "XMIPP")
                {
                    ref = initReference(data);
                    imageReal<int> mask = circularMask(ref.shape, ref.shape[0] / 2);
                    ref = ref * mask;
                }
                else if (params["-r"] == "RFAA")
                {
                    ref = reference(data, alignFFTX);
                }
                else
                {
                    string pathRef = params["-r"];
                    mrcFile mrcRef = mrcFile();
                    mrcRef.read(pathIn);
                    ref = mrcRef.data.pieceGet(0);
                }

                stackReal<double> stkBackup;
                if (XMIPP)
                    stkBackup = data;

                for (int i = 0; i < 2; ++i)
                    refinement(data, ref, alignFFTX, updateFlag);

                ref = data.stackAvg();

                if (XMIPP)
                    outInfo(stkBackup, ref, alignFFTX, info, name);
            }
            else
                throw baseException("Error: Unsupported method type! Use -h to see commandline instructions.");

            mrcFile mrc_new = mrcFile();
            stackReal<double> stk = image2Stack(ref);
            mrc_new.setData(stk);
            mrc_new.write(name + "_ref.mrcs");
            mrc_new.setData(data);
            mrc_new.write(pathOut);
            et = clock();
            std::cout << "Align time:" << (double)(et - st) / CLOCKS_PER_SEC << std::endl;
        }
    }

    return 0;
}
