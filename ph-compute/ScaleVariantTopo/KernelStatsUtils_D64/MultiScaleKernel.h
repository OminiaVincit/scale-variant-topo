#pragma once
#include "stdafx.h"
#include "KernelStatsUtils_D64.h"
#include "PersistenceUtils_D64/PersistenceBarcodes.h"

using namespace NPersistenceUtils;

namespace NMultiScaleKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const SHolePtrTypeVector<Dtype>& barcodes_1,
        const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype time);

    template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float time);

    template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(double& result, const SHolePtrTypeVector<double>& barcodes_1,
        const SHolePtrTypeVector<double>& barcodes_2, double time);

	template <class Dtype>
	KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const SHolePtrTypeVector<Dtype>& barcodes_1,
        const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype T1, Dtype T2);

	template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float T1, float T2);

	template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(double& result, const SHolePtrTypeVector<double>& barcodes_1,
        const SHolePtrTypeVector<double>& barcodes_2, double T1, double T2);

    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(Dtype& result, const SHolePtrTypeVector<Dtype>& barcodes_1,
        const SHolePtrTypeVector<Dtype>& barcodes_2, Dtype time);

    template KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float time);

    template KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(double& result, const SHolePtrTypeVector<double>& barcodes_1,
        const SHolePtrTypeVector<double>& barcodes_2, double time);

	KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float time);

    KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const SHolePtrTypeVector<float>& barcodes_1,
        const SHolePtrTypeVector<float>& barcodes_2, float T1, float T2);
}