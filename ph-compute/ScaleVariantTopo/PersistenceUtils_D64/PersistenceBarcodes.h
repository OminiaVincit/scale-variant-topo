#pragma once
#include "PersistenceUtils_D64.h"
#include <vector>

namespace NPersistenceUtils {
    template <typename Dtype>
    struct PersistenceUtils_D64_API SHole {
        SHole() {};
        SHole(Dtype _birth, Dtype _death, Dtype _tau) :
            birth(_birth), death(_death), tau(_tau) {};
        Dtype birth = 0;
        Dtype death = 0;
        Dtype tau   = 0;
        Dtype extra = 0;
    };

    template <typename Dtype>
    using SHolePtrType = std::shared_ptr<SHole<Dtype>>;

    template <typename Dtype>
    using SHolePtrTypeVector = std::vector<SHolePtrType<Dtype>>;

    template <typename Dtype>
    struct PersistenceUtils_D64_API SHoleParam {
        SHoleParam() {};
        SHoleParam(unsigned int _dim, Dtype _valInf, bool _skipInf, Dtype _thres) :
            dim(_dim), valInf(_valInf), skipInf(_skipInf), threshold(_thres) {};
        unsigned int dim = 0;
        Dtype valInf = 0;
        bool skipInf = true;
        Dtype threshold = 0;
    };

    template <typename Dtype>
    struct PersistenceUtils_D64_API SKerParam {
        SKerParam() {};
        SKerParam(Dtype _singletau, Dtype _maxtau, Dtype _interval) :
            single_tau(_singletau), max_tau(_maxtau), interval(_interval) {};
        Dtype single_tau = 0;
        Dtype max_tau = 0;
		Dtype interval = 1;
    };

    template <typename Dtype>
    class PersistenceUtils_D64_API CPersistenceBarcodes {
    public:
        CPersistenceBarcodes() {};
        CPersistenceBarcodes(const std::wstring fileName, const SHoleParam<Dtype> rprms, const SKerParam<Dtype> kprms);
        CPersistenceBarcodes(const std::wstring fileName, const SHoleParam<Dtype> rprms);
        ~CPersistenceBarcodes();
        size_t dim() const { return m_dim; };
        std::wstring GetFile() { return m_filesrc; };

        SHolePtrTypeVector<Dtype> barcodes() { return m_barcodes; };
        const size_t numbars() { return m_barcodes.size(); };
        const bool IsEmpty() { return m_barcodes.empty(); };

        std::vector<Dtype> births();
        std::vector<Dtype> deaths();

        bool PushBar(const Dtype birth, const Dtype death) {
            SHolePtrType<Dtype> hole_ptr(new SHole<Dtype>(birth, death, (Dtype)0.0));
            m_barcodes.push_back(hole_ptr);
            return true;
        }

        bool PushBar(const Dtype birth, const Dtype death, const Dtype tau) {
            SHolePtrType<Dtype> hole_ptr(new SHole<Dtype>(birth, death, tau));
            m_barcodes.push_back(hole_ptr);
            return true;
        }

        bool SetDim(size_t _dim) {
            m_dim = _dim;
            return true;
        }

        std::wstring ToFile(const std::wstring outpath, const std::wstring fileInitial);
        Dtype GetOptimalTimeHole();
    private:
        SHolePtrTypeVector<Dtype> m_barcodes = {};
        size_t m_dim = 0;
        std::wstring m_filesrc = L"";
    };
    
    template <typename Dtype>
    using CPersistenceBarcodesPtr = std::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>>;

    template <typename Dtype>
    using CPersistenceBarcodesCPtr = std::shared_ptr<const NPersistenceUtils::CPersistenceBarcodes<Dtype>>;

    template <typename Dtype>
    using TypeBarcodesPtrVec = std::vector<CPersistenceBarcodesPtr<Dtype>>;
}

namespace NPersistenceUtils {
    PersistenceUtils_D64_API bool dummy_func();

    template <typename Dtype>
    PersistenceUtils_D64_API bool MakeDiagramVec(TypeBarcodesPtrVec<Dtype> &diagram_vec,
        const std::wstring barcodes_path,
        const SHoleParam<Dtype> rprms,
        const SKerParam<Dtype>  kprms);

    template PersistenceUtils_D64_API bool MakeDiagramVec(TypeBarcodesPtrVec<float> &diagram_vec,
        const std::wstring barcodes_path,
        const SHoleParam<float> rprms,
        const SKerParam<float>  kprms);

    template PersistenceUtils_D64_API bool MakeDiagramVec(TypeBarcodesPtrVec<double> &diagram_vec,
        const std::wstring barcodes_path,
        const SHoleParam<double> rprms,
        const SKerParam<double>  kprms);
}