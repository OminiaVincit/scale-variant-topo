// -------------------------------------------------------------------
// $ProjectName     : DiagramDistanceRunner_D64
// $FileName        : DiagramDistanceRunner_D64.cpp$
// $Programmer      : Tran Quoc Hoan$
// -------------------------------------------------------------------

#include "stdafx.h"
#include "fstream"
#include <ostream>
#include <filesystem>

#include "DiagramDistanceRunner_D64.h"
#include "../KernelStatsUtils_D64/KernelFisherDA.h"
#include "../KernelStatsUtils_D64/MultiScaleKernel.h"
#include "../KernelStatsUtils_D64/RiemannianManifoldKernel.h"
#include "../PersistenceUtils_D64/PersistenceBarcodes.h"
#include "../PersistenceUtils_D64/PersistenceDeclarations.h"
#include "../PersistenceUtils_D64/PersistenceUtils_D64.h"

#include "../TopoUtils_D64/FileUtils.h"
#include "../TopoUtils_D64/StringUtils.h"

#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <Eigen/Dense>

using namespace boost::program_options;
using namespace NPersistenceUtils;
using namespace NKernelStatsUtils;
using namespace Eigen;
namespace fs = std::experimental::filesystem;

namespace NDiagramDistanceRunner {
    using namespace NMultiScaleKernelUtils;
    typedef Eigen::Matrix<FType, Dynamic, Dynamic> MaxtrixFType;

    enum KernelType {
        L2_INNER_MULTI_SCALE_SSE,
        L2_SQUARE_DIST,
        SLICE_WASS,
        L2_INNER_MULTI_SCALE_NOSSE,
        RIEMANNIAN_METRIC,
    };

    template <class Dtype>
    struct KernelPrm {
        KernelType kertype = KernelType::L2_INNER_MULTI_SCALE_SSE;
        Dtype T1 = (Dtype)1.0;
        Dtype T2 = (Dtype)1.0; // T_2 = time_hole * time_rate
        size_t theta_ndirs = 1;
        size_t phi_ndirs = 1;
        Dtype gamma = (Dtype)1.0;
        std::wstring postfix = L"";
    };

    bool KernelProduct(FType &result, CPersistenceBarcodesPtr<FType> &bar1, CPersistenceBarcodesPtr<FType> &bar2, KernelPrm<FType>& prm) {
        auto T1 = prm.T1;
        auto T2 = prm.T2;

        switch (prm.kertype)
        {
        case KernelType::L2_INNER_MULTI_SCALE_SSE:
            L2DiagramMskInnerProductSSE(result, bar1->barcodes(), bar2->barcodes(), T1, T2);
            break;
        case KernelType::L2_SQUARE_DIST:
            L2DiagramMskDistanceSquare(result, bar1->barcodes(), bar2->barcodes(), T1);
            break;
        case KernelType::L2_INNER_MULTI_SCALE_NOSSE:
            L2DiagramMskInnerProduct(result, bar1->barcodes(), bar2->barcodes(), T1, T2);
            break;
        case KernelType::RIEMANNIAN_METRIC:
            NRiemannianManifoldKernelUtils::RiemannGeodesicMetric(result, bar1->barcodes(), bar2->barcodes(), T1, T2);
            break;
        default:
            break;
        }
        return true;
    }

    MaxtrixFType KernelVecsProduct(TypeBarcodesPtrVec<FType> diagram_vecs, KernelPrm<FType>& prm) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        size_t num_diagrams = diagram_vecs.size();
        MaxtrixFType gram_mat(num_diagrams, num_diagrams);

        for (size_t i = 0; i < num_diagrams; ++i) {
            auto total = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - begin).count();
            auto hours = total / 3600;
            auto mins = (total - hours * 3600) / 60;
            auto secs = total - hours * 3600 - mins * 60;
            std::cout << "Ellapsed time (h:m:s)= " << hours << ":" << mins << ":" << secs << 
                ", Index of diagram: " << i << " per total of " << num_diagrams << std::endl;

            // multi-threads
            concurrency::parallel_for(i, num_diagrams, [&](size_t j) {
                FType result(0);
                KernelProduct(result, diagram_vecs[i], diagram_vecs[j], prm);
                gram_mat(i, j) = result;
                //std::cout << i << " " << j << " " << diagram_vecs[i]->numbars() 
                //    << " " << diagram_vecs[j]->numbars() << std::endl;
                return true;
            });
        }
        // gram matrix should be symmetric
        for (int i = 0; i < num_diagrams; ++i) {
            for (int j = 0; j < i; ++j) {
                gram_mat(i, j) = gram_mat(j, i);
            }
        }
        return gram_mat;
    }

    MaxtrixFType KernelVecsProduct(TypeBarcodesPtrVec<FType> dvecs1, TypeBarcodesPtrVec<FType> dvecs2, KernelPrm<FType>& prm) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        size_t nums1 = dvecs1.size();
        size_t nums2 = dvecs2.size();
        if (nums1 == 0) return KernelVecsProduct(dvecs2, prm);
        if (nums2 == 0) return KernelVecsProduct(dvecs1, prm);

        MaxtrixFType gram_mat(nums1, nums2);

        for (size_t i = 0; i < nums1; ++i) {
            auto total = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - begin).count();
            auto hours = total / 3600;
            auto mins = (total - hours * 3600) / 60;
            auto secs = total - hours * 3600 - mins * 60;
            std::cout << "Ellapsed time (h:m:s)= " << hours << ":" << mins << ":" << secs <<
                ", Index of diagram: " << i << " per total of " << nums1 << std::endl;
            
            // multi-threads
            concurrency::parallel_for((size_t)0, nums2, [&](size_t j) {
                FType result(0);
                KernelProduct(result, dvecs1[i], dvecs2[j], prm);
                gram_mat(i, j) = result;
                return true;
            });
        }
        return gram_mat;
    }

    // Get the median of an unordered set of numbers of arbitrary 
    // type without modifying the underlying dataset.
    template <typename It>
    auto GetMedian(It begin, It end)
    {
        using T = typename std::iterator_traits<It>::value_type;
        std::vector<T> data(begin, end);
        std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
        return data[data.size() / 2];
    }

    bool FindOptimalTimeHole(TypeBarcodesPtrVec<FType> dvecs, FType &tHole) {
        if(dvecs.empty())
            return false;
        auto N = dvecs.size();
        std::vector<FType> time_hole(N, (FType)0.0);
        concurrency::parallel_for((size_t)0, N, [&](size_t z) {
            time_hole[z] = dvecs[z]->GetOptimalTimeHole();
        });

        tHole = GetMedian(time_hole.begin(), time_hole.end());
        return true;
    }
    bool FindOptimalTimeHole(TypeBarcodesPtrVec<FType> dvecs1, TypeBarcodesPtrVec<FType> dvecs2, FType &tHole) {
        TypeBarcodesPtrVec<FType> dvecs = dvecs1;
        dvecs.insert(dvecs.end(), dvecs2.begin(), dvecs2.end());
        FindOptimalTimeHole(dvecs, tHole);
        return true;
    }

    bool MakeGramMatAndSaveResultToFile(TypeBarcodesPtrVec<FType> dvecs1, TypeBarcodesPtrVec<FType> dvecs2,
        std::wstring& output_path, KernelPrm<FType>& prm) {
        std::cout << "Computing gram matrix and save to file: " << NStringUtil::_w2s(output_path) << std::endl;
        fs::path opath(output_path);
        if (prm.kertype == KernelType::L2_INNER_MULTI_SCALE_NOSSE ||
            prm.kertype == KernelType::L2_INNER_MULTI_SCALE_SSE ||
            prm.kertype == KernelType::L2_SQUARE_DIST ||
            prm.kertype == KernelType::RIEMANNIAN_METRIC) {
            if (prm.T1 == 0.0) FindOptimalTimeHole(dvecs1, dvecs2, prm.T1);
            std::cout << "Optimal T1 = " << prm.T1 << ", T2 = " << prm.T2 << std::endl;
        }

        MaxtrixFType gram_mat = KernelVecsProduct(dvecs1, dvecs2, prm);

        const size_t rows = gram_mat.rows();
        const size_t cols = gram_mat.cols();

        std::wstring outfile = output_path;
        if (fs::is_directory(opath)) {
            fs::path file(prm.postfix + L".txt");
            outfile = (opath / file).c_str();
        }

        // create folder if not exist
        fs::path p(outfile);
        fs::path dir = p.parent_path();
        auto dst_folder = dir.wstring();
        if (!fs::exists(dst_folder))
            fs::create_directory(dst_folder);

        // write result
        std::ofstream outf(outfile);
        Eigen::IOFormat fmt(FullPrecision, DontAlignCols, "\t");
        if (outf.is_open()) {
            outf << gram_mat.format(fmt) << std::endl;
        }
        outf.close();
        return true;
    }

    bool SaveKernelResultFromBarcodeListFile(std::wstring& barcodes_path, std::wstring& output_path, KernelPrm<FType>& prm,
        const SHoleParam<FType> rprms, const SKerParam<FType> kpm) {
        TypeBarcodesPtrVec<FType> dvecs1, dvecs2;
        MakeDiagramVec(dvecs1, barcodes_path, rprms, kpm);
        if(dvecs1.empty()) 
            return false;
        MakeGramMatAndSaveResultToFile(dvecs1, dvecs2, output_path, prm);
        return true;
    }

    bool SaveKernelResultFromBarcodeListFile(std::wstring& barcodes_lpath, std::wstring& barcodes_rpath,
        std::wstring& output_path, KernelPrm<FType>& prm,
        const SHoleParam<FType> rprms, const SKerParam<FType> kpm) {
        TypeBarcodesPtrVec<FType> dvecs1, dvecs2;
        std::cout << "Making barcodes for left path" << std::endl;
        MakeDiagramVec(dvecs1, barcodes_lpath, rprms, kpm);
        std::cout << "Making barcodes for right path" << std::endl;
        MakeDiagramVec(dvecs2, barcodes_rpath, rprms, kpm);
        if(dvecs1.empty() == true && dvecs2.empty() == true)
            return false;
        MakeGramMatAndSaveResultToFile(dvecs1, dvecs2, output_path, prm);
        return true;
    }

    bool SaveKFDRFromGrammatToFile(const std::wstring& outfile, const MaxtrixFType gram_mat, std::vector<double> gammas) {
        // create folder if not exist
        fs::path p(outfile);
        fs::path dir = p.parent_path();
        auto dst_folder = dir.wstring();
        if (!fs::exists(dst_folder))
            fs::create_directory(dst_folder);

        namespace io = boost::iostreams;
        io::stream_buffer<io::file_sink> buf(NStringUtil::_w2s(outfile));
        std::ostream out(&buf);
        for (auto &gamma : gammas) {
            auto kfdrs = ComputeKFDRs((int)gram_mat.rows(), gram_mat, gamma);
            out << gamma << ' ';
            for (auto val : kfdrs) {
                out << val << ' ';
            }
            out << std::endl;
        }
        return true;
    }


    bool SaveKFDRFromBarcodeListFile(std::wstring& barcodes_path, std::wstring& kernel_outfile, KernelPrm<FType>& prm,
        const SHoleParam<FType>& rprms, const SKerParam<FType>& kpm, std::vector<double> gammas) {
        TypeBarcodesPtrVec<FType> dvecs1, dvecs2;
        MakeDiagramVec(dvecs1, barcodes_path, rprms, kpm);
        if(dvecs1.empty())
            return false;

        if (prm.kertype == KernelType::L2_INNER_MULTI_SCALE_NOSSE ||
            prm.kertype == KernelType::L2_INNER_MULTI_SCALE_SSE ||
            prm.kertype == KernelType::L2_SQUARE_DIST ||
            prm.kertype == KernelType::RIEMANNIAN_METRIC) {
            if (prm.T1 == 0.0) FindOptimalTimeHole(dvecs1, dvecs2, prm.T1);
        }

        MaxtrixFType gram_mat = KernelVecsProduct(dvecs1, dvecs2, prm);
        // save gram_mat file
        // write result
        std::ofstream outf(kernel_outfile);
        Eigen::IOFormat fmt(FullPrecision, DontAlignCols, "\t");
        if (outf.is_open()) {
            outf << gram_mat.format(fmt) << std::endl;
        }
        outf.close();
        fs::path p(kernel_outfile);
        std::wstring kfdr_outfile = p.parent_path().c_str();
        kfdr_outfile = kfdr_outfile + L"/kfdr_" + p.filename().c_str(); // des.wstring();
        
        return SaveKFDRFromGrammatToFile(kfdr_outfile, gram_mat, gammas);
    }

    bool ReadKernelFileFromFileList(std::vector<std::wstring>&kernel_list, std::wstring& kernel_list_file)
    {
        std::wifstream read_op;
        read_op.open(kernel_list_file);
        if (!read_op.good()) {
            return false;
        }
        else {
            while (read_op.good()) {
                std::wstring line;
                std::getline(read_op, line);

                if (!line.empty()) {
                    kernel_list.push_back(line);
                }
            }
        }
        return true;
    }
    bool SaveKFDRFromKernelList(std::wstring& kernel_list_file, std::wstring& dst_folder, size_t nbegin, size_t nend, std::vector<double> gammas) {
        if (!fs::exists(dst_folder))
            fs::create_directory(dst_folder);

        std::vector<std::wstring> kernel_list;
        ReadKernelFileFromFileList(kernel_list, kernel_list_file);
        for (auto kf : kernel_list) {
            std::vector<FType> xs;
            std::wifstream in(kf);
            std::wstring line;
            if (in.is_open()) {
                while (std::getline(in, line)) {
                    std::wstringstream wstream(line);
                    while (!wstream.eof()) {
                        FType tmp = (FType)0.0;
                        wstream >> tmp;
                        xs.push_back(tmp);
                    }
                }
            }
            in.close();
            size_t rows = (size_t)(std::sqrt(xs.size()));
            size_t cols = (size_t)(std::sqrt(xs.size()));
            if (rows*cols != xs.size())
                continue;
            size_t rs = rows;
            size_t cs = cols;
            if (nbegin < nend) {
                rs = std::min(rows, nend - nbegin + 1);
                cs = rs;
            }
            MaxtrixFType grammat = MaxtrixFType(rs, cs);
            for (size_t i = 0; i < rows; ++i) {
                if (i < nbegin) continue;
                if (i >= nbegin + rs) continue;
                for (size_t j = 0; j < cols; ++j) {
                    if (j < nbegin) continue;
                    if (j >= nbegin + rs) continue;
                    grammat(i - nbegin, j - nbegin) = xs[i*cols + j];
                }
            }
            fs::path p(kf);
            std::wstring outfile = dst_folder + L"/kfdr_" + p.filename().c_str();
            if (nbegin < nend)
                outfile = dst_folder + L"/kfdr_bg_" + std::to_wstring(nbegin) + L"_ed_"+ std::to_wstring(nend) + L"_" + p.filename().c_str(); // des.wstring();
            return SaveKFDRFromGrammatToFile(outfile, grammat, gammas);
        }
        return true;
    }
}

int main(int argc, char** argv)
{
    using namespace NDiagramDistanceRunner;
    int nRetCode = 0;

    options_description description("DiagramDistance");
    description.add_options()
        ("T1", value<double>()->default_value(0.0), "Parameter T1 in the kernel (=0 for optimal value)")
        ("T2", value<double>()->default_value(1.0),  "Parameter T2 in the kernel")
        ("left,l", value<std::string>()->default_value(""), "(left) Input as List of barcodes")
        ("right,r", value<std::string>()->default_value(""), "(right) Input as List of barcodes")
        ("dim,d", value<unsigned>()->default_value(0), "Dimension of holes to compute kernel")
        ("infval", value<double>()->default_value(0.0), 
            "If infval<=0.0 skip holes with infinity death-scale, otherwise replace these death-scales with a default value")
        ("thres", value<double>()->default_value(0.0), "Threshold to skip holes with death-birth < thres (default=0 to use all holes)")
        ("output,o", value<std::string>()->default_value("grammat.txt"), "Output file of gram matrix")
        ("posfix", value<std::string>()->default_value(""), "postfix for output file")
        ("method", value<int>()->default_value(0),
            "0: L2_inner_multi_sse, 1: L2_squared_distance, 2: Slice Wasserstein, 3:L2_inner_multi_nosse, 4:riemmannian_metric")
        ("single_tau", value<double>()->default_value(-1.0), "Specify for single tau, < 0: take all tau")
        ("max_tau", value<double>()->default_value(0.0), 
            "Maximum of scale when calculating kernel, max_tau = 0 means that taking all possible tau")
		("interval", value<double>()->default_value(1.0),
				"Interval to skip tau when calculating kernel with all taus")
        ("kfdr", value<bool>()->default_value(false), "option to calculate kfdr, 1: kfdr, 0:kernel")
        ("kfdrout", value<std::string>()->default_value("kfdr"), "Output folder for kfdr")
        ("nbegin", value<unsigned>()->default_value(0))
        ("nend", value<unsigned>()->default_value(0))
        ("kerls", value<std::string>()->default_value(""), "Input as list of kernel for kfdr")
        ("help,H", "Help: Diagram Distance to compute kernel of diagrams")
        ("version,v", "v1.0")
        ;
    variables_map vm;
    store(parse_command_line(argc, argv, description), vm);
    notify(vm);
    if (vm.count("help")) {
        std::cout << description << std::endl;
        return nRetCode;
    }
    
    auto T1 = static_cast<FType>(vm["T1"].as<double>());
    auto T2 = static_cast<FType>(vm["T2"].as<double>());

    auto nbegin = static_cast<size_t>(vm["nbegin"].as<unsigned>());
    auto nend = static_cast<size_t>(vm["nend"].as<unsigned>());
    auto dim = vm["dim"].as<unsigned>();
    auto single_tau = vm["single_tau"].as<double>();
    auto max_tau = static_cast<FType>(vm["max_tau"].as<double>());
	auto interval = static_cast<FType>(vm["interval"].as<double>());
	if (interval <= 0.0) interval = 1.0;

    bool skip = true;
    auto infval = static_cast<FType>(vm["infval"].as<double>());
    if (infval > 0.0) skip = false;

    auto threshold = static_cast<FType>(vm["thres"].as<double>());
    
    auto kernel_list_path = NStringUtil::_s2w(vm["kerls"].as<std::string>());
    auto left_path = NStringUtil::_s2w(vm["left"].as<std::string>());
    auto right_path = NStringUtil::_s2w(vm["right"].as<std::string>());

    auto output_kfdr = NStringUtil::_s2w(vm["kfdrout"].as<std::string>());
    auto output_path = NStringUtil::_s2w(vm["output"].as<std::string>());
    auto posfix = NStringUtil::_s2w(vm["posfix"].as<std::string>());

    auto method = vm["method"].as<int>();

    KernelType kertype;
    switch (method) {
    case 0:
        kertype = KernelType::L2_INNER_MULTI_SCALE_SSE;
        break;
    case 1:
        kertype = KernelType::L2_SQUARE_DIST;
        break;
    case 2:
        kertype = KernelType::SLICE_WASS;
        break;
    case 3:
        kertype = KernelType::L2_INNER_MULTI_SCALE_NOSSE;
        break;
    case 4:
        kertype = KernelType::RIEMANNIAN_METRIC;
        break;
    default:
        std::cout << "Kernel method is not defined!" << std::endl;
        std::cout << description;
        return nRetCode;
    }
    
    KernelPrm<FType> prm;
    prm.kertype = kertype;
    prm.T1 = T1;
    prm.T2 = T2 * interval * interval;
    prm.postfix = posfix + L"_T1_" + std::to_wstring(T1) +
        L"_T2_" + std::to_wstring(T2);

    const SHoleParam<FType> rprms(dim, infval, skip, threshold);
    const SKerParam<FType> kpm(single_tau, max_tau, interval);

    auto kfdr = vm["kfdr"].as<bool>();
    if (kfdr == TRUE) {
        std::vector<double> gammas = {
            1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6,
            1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5,
            1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
            1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
            1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2,
            1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1,
            1e+0, 2e+0, 3e+0, 4e+0, 5e+0, 6e+0, 7e+0, 8e+0, 9e+0,
            1e+1, 2e+1, 3e+1, 4e+1, 5e+1, 6e+1, 7e+1, 8e+1, 9e+1,
            1e+2, 2e+2, 3e+2, 4e+2, 5e+2, 6e+2, 7e+2, 8e+2, 9e+2,
            1e+3, 2e+3, 3e+3, 4e+3, 5e+3, 6e+3, 7e+3, 8e+3, 9e+3,
            1e+4, 2e+4, 3e+4, 4e+4, 5e+4, 6e+4, 7e+4, 8e+4, 9e+4,
            1e+5, 2e+5, 3e+5, 4e+5, 5e+5, 6e+5, 7e+5, 8e+5, 9e+5
        };
        if (kernel_list_path == L"") {
            SaveKFDRFromBarcodeListFile(left_path, output_path, prm, rprms, kpm, gammas);
        }
        else {
            SaveKFDRFromKernelList(kernel_list_path, output_kfdr, nbegin, nend, gammas);
        }
    }
    else {
        if (left_path != L"" && !fs::exists(left_path)) {
            std::cout << "Not found left path for calculating kernel: " << NStringUtil::_w2s(left_path) << std::endl;
            return nRetCode;
        }
        if (right_path != L"" && !fs::exists(right_path)) {
            std::cout << "Not found right path for calculating kernel: " << NStringUtil::_w2s(right_path) << std::endl;
            return nRetCode;
        }
        SaveKernelResultFromBarcodeListFile(left_path, right_path, output_path, prm, rprms, kpm);
    }
    return nRetCode;
}
