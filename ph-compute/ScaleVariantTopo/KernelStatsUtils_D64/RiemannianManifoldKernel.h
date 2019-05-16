#pragma once

#include "stdafx.h"
#include "KernelStatsUtils_D64.h"
#include "PersistenceUtils_D64/PersistenceBarcodes.h"

using namespace NPersistenceUtils;

namespace NRiemannianManifoldKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool RiemannGeodesicMetric(Dtype& result, const SHolePtrTypeVector<Dtype>& barcodes_1,
        const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype T1, Dtype T2);

    KernelStatsUtils_D64_API bool RiemannGeodesicMetric(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float T1, float T2);
}