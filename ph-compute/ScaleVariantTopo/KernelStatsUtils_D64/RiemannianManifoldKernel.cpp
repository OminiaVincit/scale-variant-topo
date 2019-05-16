#include "stdafx.h"
#include "RiemannianManifoldKernel.h"

#include <numeric>
#include <ppl.h>
#include <ppltasks.h>

static const double THRES_ZERO = 70;

template <class Dtype>
bool ProjectToDiagonalAndAppend(SHolePtrTypeVector<Dtype>& barcodes_1,
    SHolePtrTypeVector<Dtype>& barcodes_2) {
    SHolePtrTypeVector<Dtype> proj_1to2;
    for (auto bar : barcodes_1) {
        auto tmp = (bar->birth + bar->death) / 2.0;
        SHolePtrType<Dtype> barptr(new SHole<Dtype>(tmp, tmp, bar->tau));
        proj_1to2.push_back(barptr);
    }

    for (auto bar : barcodes_2) {
        auto tmp = (bar->birth + bar->death) / 2.0;
        SHolePtrType<Dtype> barptr(new SHole<Dtype>(tmp, tmp, bar->tau));
        barcodes_1.push_back(barptr);
    }
    std::copy(proj_1to2.begin(), proj_1to2.end(), std::back_inserter(barcodes_2));
    return true;
}

template <class Dtype>
Dtype GetNorm(SHolePtrType<Dtype>& bar1, SHolePtrType<Dtype>& bar2, const Dtype& time) {
    auto diff_1 = bar1->birth - bar2->birth;
    auto diff_2 = bar1->death - bar2->death;
    auto diff_3 = bar1->tau   - bar2->tau;
    auto diff_4 = bar1->extra - bar2->extra;
    auto diff_square = diff_1*diff_1 + diff_2*diff_2 + diff_3*diff_3 + diff_4 * diff_4;
    const Dtype result = diff_square / (2.0*time);
    return result;
}

float GetNormSSE(SHolePtrType<float>& bar1, SHolePtrType<float>& bar2, const float& time) {
    float result = (float)0.0;
    __m128 diff = _mm_sub_ps(
        _mm_set_ps(bar1->extra, bar1->tau, bar1->death, bar1->birth),
        _mm_set_ps(bar2->extra, bar2->tau, bar2->death, bar2->birth)
    );
    diff = _mm_mul_ps(diff, diff);
    diff = _mm_hadd_ps(diff, diff);
    diff = _mm_hadd_ps(diff, diff);
    const float val = diff.m128_f32[0] / (2.0*time);
    return val;
}
template <class Dtype>
bool RiemannGeodesicMetricSingle(Dtype& result, const SHolePtrTypeVector<Dtype>& barcodes_1,
    const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype time)
{
    auto tbars_1 = barcodes_1;
    auto tbars_2 = barcodes_2;

    ProjectToDiagonalAndAppend(tbars_1, tbars_2);
    if(tbars_1.size() != tbars_2.size())
        return false;

    auto N = tbars_1.size();
    auto m = tbars_1.size() + tbars_2.size();
    std::vector<Dtype> vx(m, (Dtype)0.0);
    std::vector<Dtype> vy(m, (Dtype)0.0);
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N; ++i) {
            auto diff_11 = GetNorm<Dtype>(tbars_1[i], tbars_1[j], time);
            auto diff_12 = GetNorm<Dtype>(tbars_1[i], tbars_2[j], time);
            auto diff_21 = GetNorm<Dtype>(tbars_2[i], tbars_1[j], time);
            auto diff_22 = GetNorm<Dtype>(tbars_2[i], tbars_2[j], time);

            if (diff_11 < THRES_ZERO) vx[j]     += std::exp(-diff_11);
            if (diff_12 < THRES_ZERO) vx[j + N] += std::exp(-diff_12);
            if (diff_21 < THRES_ZERO) vy[j]     += std::exp(-diff_21);
            if (diff_22 < THRES_ZERO) vy[j + N] += std::exp(-diff_22);
        }
    }
    const Dtype x_sum = std::accumulate(vx.begin(), vx.end(), 0.0, std::plus<Dtype>());
    const Dtype y_sum = std::accumulate(vy.begin(), vy.end(), 0.0, std::plus<Dtype>());
    for (auto &x : vx) x = std::sqrt(x / x_sum);
    for (auto &y : vy) y = std::sqrt(y / y_sum);

    Dtype product = 0;
    for (size_t i = 0; i < m; ++i) product += vx[i] * vy[i];
    product = std::max<Dtype>(product, 0.0);
    product = std::min<Dtype>(product, 1.0);
    
    auto val = std::acos(product);
    result = val;
    return true;
}

template <class Dtype>
KernelStatsUtils_D64_API bool NRiemannianManifoldKernelUtils::RiemannGeodesicMetric(Dtype & result, 
    const SHolePtrTypeVector<Dtype>& barcodes_1,
    const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype T1, Dtype T2)
{
    if (T1 <= 0 || T2 <= 0)
        return false;
    const Dtype rate = std::sqrt(T1 / T2);
    SHolePtrTypeVector<Dtype> barvec1;
    SHolePtrTypeVector<Dtype> barvec2;
    for (auto tmp : barcodes_1) {
        SHolePtrType<Dtype> barptr(new SHole<Dtype>(tmp->birth, tmp->death, tmp->tau));
        barvec1.push_back(barptr);
    }

    for (auto tmp : barcodes_2) {
        SHolePtrType<Dtype> barptr(new SHole<Dtype>(tmp->birth, tmp->death, tmp->tau));
        barvec2.push_back(barptr);
    }
    RiemannGeodesicMetricSingle(result, barvec1, barvec2, rate);
    return true;
}


bool RiemannGeodesicMetricSSE(float& result, 
    const SHolePtrTypeVector<float>& barcodes_1,
    const SHolePtrTypeVector<float>& barcodes_2, float time)
{
    auto tbars_1 = barcodes_1;
    auto tbars_2 = barcodes_2;

    ProjectToDiagonalAndAppend(tbars_1, tbars_2);
    if (tbars_1.size() != tbars_2.size())
        return false;

    auto N = tbars_1.size();
    auto m = tbars_1.size() + tbars_2.size();
    std::vector<float> vx(m, (float)0.0);
    std::vector<float> vy(m, (float)0.0);
    __m128 zero = _mm_set_ps1((float)0.0);

    std::vector<__m128> vss(N, _mm_set_ps1((float)0.0));
    for (size_t j = 0; j < N; ++j) {
        __m128 tmp = zero;
        for (size_t i = 0; i < N; ++i) {
            auto diff_11 = GetNormSSE(tbars_1[i], tbars_1[j], time);
            auto diff_12 = GetNormSSE(tbars_1[i], tbars_2[j], time);
            auto diff_21 = GetNormSSE(tbars_2[i], tbars_1[j], time);
            auto diff_22 = GetNormSSE(tbars_2[i], tbars_2[j], time);
            tmp = _mm_add_ps(tmp, _mm_set_ps(std::exp(-diff_22), std::exp(-diff_21), std::exp(-diff_12), std::exp(-diff_11)));
            //if (diff_11 < THRES_ZERO) vx[j] += std::exp(-diff_11);
            //if (diff_12 < THRES_ZERO) vx[j + N] += std::exp(-diff_12);
            //if (diff_21 < THRES_ZERO) vy[j] += std::exp(-diff_21);
            //if (diff_22 < THRES_ZERO) vy[j + N] += std::exp(-diff_22);
        }
        vx[j]     = tmp.m128_f32[0];
        vx[j + N] = tmp.m128_f32[1];
        vy[j]     = tmp.m128_f32[2];
        vy[j + N] = tmp.m128_f32[3];
    }

    const float x_sum = std::accumulate(vx.begin(), vx.end(), 0.0, std::plus<float>());
    const float y_sum = std::accumulate(vy.begin(), vy.end(), 0.0, std::plus<float>());
    for (auto &x : vx) x = std::sqrt(x / x_sum);
    for (auto &y : vy) y = std::sqrt(y / y_sum);
    float product = 0;
    for (size_t i = 0; i < m; ++i) product += vx[i] * vy[i];
    product = std::max<float>(product, 0.0);
    product = std::min<float>(product, 1.0);

    auto val = std::acos(product);
    result = val;
    return true;
}

bool RiemannGeodesicMetricSSE(float & result,
    const SHolePtrTypeVector<float>& barcodes_1,
    const SHolePtrTypeVector<float>& barcodes_2, float T1, float T2)
{
    if (T1 <= 0 || T2 <= 0)
        return false;
    const float rate = std::sqrt(T1 / T2);
    SHolePtrTypeVector<float> barvec1;
    SHolePtrTypeVector<float> barvec2;
    for (auto tmp : barcodes_1) {
        SHolePtrType<float> barptr(new SHole<float>(tmp->birth, tmp->death, rate * tmp->tau));
        barvec1.push_back(barptr);
    }

    for (auto tmp : barcodes_2) {
        SHolePtrType<float> barptr(new SHole<float>(tmp->birth, tmp->death, rate * tmp->tau));
        barvec2.push_back(barptr);
    }
    RiemannGeodesicMetricSSE(result, barvec1, barvec2, rate);
    return true;
}

KernelStatsUtils_D64_API bool NRiemannianManifoldKernelUtils::RiemannGeodesicMetric(float & result, 
    const SHolePtrTypeVector<float>& barcodes_1,
    const SHolePtrTypeVector<float>& barcodes_2, float T1, float T2)
{
    return NRiemannianManifoldKernelUtils::RiemannGeodesicMetric<float>(result, barcodes_1, barcodes_2, T1, T2);
    //return RiemannGeodesicMetricSSE(result, barcodes_1, barcodes_2, T1, T2);
};
