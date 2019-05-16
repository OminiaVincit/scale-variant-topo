#include "stdafx.h"
#include <sstream>
#include <fstream>
#include <ostream>
#include <filesystem>
#include "PersistenceBarcodes.h"
#include "../TopoUtils_D64/StringUtils.h"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

namespace NPersistenceUtils {
    namespace fs = std::experimental::filesystem;

    template <typename Dtype>
    CPersistenceBarcodes<Dtype>::~CPersistenceBarcodes() {
    };

    template <typename Dtype>
    CPersistenceBarcodes<Dtype>::CPersistenceBarcodes(const std::wstring filename,
        const SHoleParam<Dtype> rprms, const SKerParam<Dtype> kprms) {
        auto dim = rprms.dim;
        auto valInfty = rprms.valInf;
        auto skipInfty = rprms.skipInf;
        auto threshold = rprms.threshold;

        auto single_tau = kprms.single_tau;
        auto max_tau    = kprms.max_tau;

        //std::cout << "Filename = " << NStringUtil::_w2s(filename) << std::endl;
        std::wifstream read_op;
        read_op.open(NStringUtil::_w2s(filename));
        if (!read_op.good()) {
            return;
        }
        else {
            while (read_op.good()) {
                std::wstring line;
                std::getline(read_op, line);

                if (line.size()) {
                    std::wstringstream s;
                    std::wstring wbirth, wdeath, wtau;
                    s << line;
                    s >> wbirth;
                    s >> wdeath;
                    s >> wtau;

                    Dtype birth, death, tau;
                    if (typeid(Dtype) == typeid(double)) {
                        birth = _wtof(wbirth.c_str());
                        death = _wtof(wdeath.c_str());
                        tau   = _wtof(wtau.c_str());
                    }
                    else {
                        birth = (float) _wtof(wbirth.c_str());
                        death = (float) _wtof(wdeath.c_str());
                        tau   = (float) _wtof(wtau.c_str());
                    }
                    
                    // read with specified tau
                    if (max_tau > 0 && tau > max_tau)
                        break;
                    if (single_tau >= 0 && single_tau != tau)
                        continue;

                    if (skipInfty && death == std::numeric_limits<Dtype>::infinity()) {
                        continue;
                    }
                    if (death == std::numeric_limits<Dtype>::infinity()) {
                        death = std::max(birth, valInfty);
                    }
                    if (abs(death - birth) <= threshold) continue;
                    PushBar(birth, death, tau);
                }
                else {
                    break;
                }
            }
        }
        read_op.close();
        //std::cout << "Number of holes " << m_barcodes.size() << std::endl;
        m_dim = static_cast<size_t>(dim);
        m_filesrc = filename;
    }

    template <typename Dtype>
    CPersistenceBarcodes<Dtype>::CPersistenceBarcodes(const std::wstring filename,
        const SHoleParam<Dtype> rprms) {
        const SKerParam<Dtype> kprms(-1.0, 0.0);
        CPersistenceBarcodes<Dtype>(filename, rprms, kprms);
    }

    template <typename Dtype>
    std::vector<Dtype> CPersistenceBarcodes<Dtype>::births() {
        std::vector<Dtype> birth_vecs;
        for (auto bar : m_barcodes) birth_vecs.push_back(bar->birth);
        return birth_vecs;
    }

    template <typename Dtype>
    std::vector<Dtype> CPersistenceBarcodes<Dtype>::deaths() {
        std::vector<Dtype> death_vecs;
        for (auto bar : m_barcodes) death_vecs.push_back(bar->death);
        return death_vecs;
    }

    template<typename Dtype>
    std::wstring CPersistenceBarcodes<Dtype>::ToFile(const std::wstring outpath, const std::wstring fileInitial)
    {
        namespace io = boost::iostreams;

        std::wstring outfile = outpath + L"/barcode_" + fileInitial + L"_dim_" + std::to_wstring(dim()) + L".txt";
        io::stream_buffer<io::file_sink> buf(NStringUtil::_w2s(outfile));
        std::ostream out(&buf);
        for (auto bar : m_barcodes) {
            out << bar->birth << " " << bar->death << std::endl;
        }
        return outfile;
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

    template<typename Dtype>
    Dtype CPersistenceBarcodes<Dtype>::GetOptimalTimeHole()
    {
        std::vector<Dtype> dis_diffs;
        Dtype tHole = Dtype(0.0);
        auto n = m_barcodes.size();
        for (size_t i = 0; i < n; ++i) {
            auto bar1 = m_barcodes[i];
            for (size_t j = i + 1; j < n; ++j) {
                auto bar2 = m_barcodes[j];
                auto xdiff = bar1->birth - bar2->birth;
                auto ydiff = bar1->death - bar2->death;
                dis_diffs.push_back(xdiff * xdiff + ydiff * ydiff);
            }
        }
        if (dis_diffs.empty() == false) {
            tHole = GetMedian(dis_diffs.begin(), dis_diffs.end());
        }
        return tHole;
    }

    template PersistenceUtils_D64_API struct SHoleParam<float>;
    template PersistenceUtils_D64_API struct SHoleParam<double>;

    template PersistenceUtils_D64_API struct SKerParam<float>;
    template PersistenceUtils_D64_API struct SKerParam<double>;

    template PersistenceUtils_D64_API class CPersistenceBarcodes<float>;
    template PersistenceUtils_D64_API class CPersistenceBarcodes<double>;
}

namespace NPersistenceUtils {
    PersistenceUtils_D64_API bool dummy_func()
    {
        return true;
    }

    template <typename Dtype>
    bool MakeDiagramVecFromBarcodesFile(TypeBarcodesPtrVec<Dtype> &diagram_vec,
        const std::wstring barlist_filename, const SHoleParam<Dtype> rprms, const SKerParam<Dtype> kprms) {
        // Read parameters from file
        std::wifstream read_op;
        read_op.open(barlist_filename);
        if (!read_op.good()) {
            return false;
        }
        else {
            while (read_op.good()) {
                std::wstring line;
                std::getline(read_op, line);

                if (!line.empty()) {
                    std::wstring bar_path = L"";
                    fs::path bpath(line);
                    if (fs::is_regular_file(bpath) == false) {
                        fs::path barlist_path(barlist_filename);
                        bpath = barlist_path.remove_filename() / bpath;
                        if (fs::is_regular_file(bpath)) bar_path = bpath.c_str();
                        else continue;
                    }
                    else {
                        bar_path = line;
                    }
                    //std::cout << "Line = " << NStringUtil::_w2s(line) << std::endl;
                    //std::cout << "barlist = " << NStringUtil::_w2s(barlist_filename) << std::endl;
                    //std::cout << "Barpath = " << NStringUtil::_w2s(bar_path) << std::endl;
                    CPersistenceBarcodesPtr<Dtype> bar_ptr(new CPersistenceBarcodes<Dtype>(bar_path, rprms, kprms));
                    diagram_vec.push_back(bar_ptr);
                }
            }
        }
        read_op.close();
        return true;
    }

    template <typename Dtype>
    PersistenceUtils_D64_API bool MakeDiagramVec(TypeBarcodesPtrVec<Dtype>& diagram_vec,
        const std::wstring barcodes_path,
        const SHoleParam<Dtype> rprms, const SKerParam<Dtype> kprms) {
        fs::path input_path(barcodes_path);

        if (fs::is_directory(input_path)) {
            // find barcode files in folder
            fs::directory_iterator it(input_path), eod;
            BOOST_FOREACH(fs::path const &p, std::make_pair(it, eod)) {
                if (fs::is_regular_file(p)) {
                    CPersistenceBarcodesPtr<Dtype> bar_ptr(new CPersistenceBarcodes<Dtype>(p.c_str(), rprms));
                    if (bar_ptr->IsEmpty() == false) {
                        diagram_vec.push_back(bar_ptr);
                    }
                }
            }
        }
        else if (fs::is_regular_file(input_path)) {
            MakeDiagramVecFromBarcodesFile(diagram_vec, barcodes_path, rprms, kprms);
        }
        else {
            return false;
        }
        return true;
    }
}