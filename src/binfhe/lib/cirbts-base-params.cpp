#include "cirbts-base-params.h"

namespace lbcrypto{
void CirBTSCryptoParams::PreCompute() {
    auto& polyParams = m_RGSWParams1->GetPolyParams();
    //Generate LUT
    auto numLUT = m_RGSWParams2->GetDigitsGA();
    auto bitwidth = static_cast<uint32_t>(std::ceil(std::log2(numLUT)));
    const uint32_t traceShift = m_RLWEParams ? m_RLWEParams->GetTraceShift() : 0;
    const uint32_t step = (traceShift == 0) ? 1u : (1u << traceShift);

    auto Q = m_RGSWParams1->GetQ();
    auto N = m_RGSWParams1->GetN();
    const auto& Gpow = m_RGSWParams2->GetAGPower();

    NativeVector LUT(N, Q, NativeInteger(0));
    const uint32_t span = bitwidth + 1 + traceShift;
    const uint32_t blocks = (span >= 31u) ? 0u : static_cast<uint32_t>(N >> span);
    for (uint32_t i = 0; i < blocks; i++){
        for (uint32_t j = 0; j < numLUT; j++){
            const uint32_t base = static_cast<uint32_t>(((i << bitwidth) + j) * step);
            LUT[base].ModSubFastEq((Gpow[j] >> 1), Q);
            LUT[(N >> 1) + base] = Gpow[j] >> 1;
        }
    }

    m_LUT = NativePoly(polyParams, Format::COEFFICIENT, true);
    m_LUT.SetValues(std::move(LUT), Format::COEFFICIENT);
    m_LUT.SetFormat(EVALUATION);

    // Computes polynomials X^{-i} that are needed in the circuit bootstrapping
    constexpr NativeInteger one{1};
    m_monomials.reserve(2 * N);
    NativePoly aPoly(polyParams, Format::COEFFICIENT, true);
    aPoly[0].ModAddFastEq(one, Q);
    aPoly.SetFormat(Format::EVALUATION);
    m_monomials.push_back(std::move(aPoly));

    for(uint32_t i = 1; i <= N; i++){
        NativePoly aPoly(polyParams, Format::COEFFICIENT, true);
        aPoly[N - i].ModSubFastEq(one, Q);
        aPoly.SetFormat(Format::EVALUATION);
        m_monomials.push_back(std::move(aPoly));
    }

    for(uint32_t i = 1; i < N; i++){
        NativePoly aPoly(polyParams, Format::COEFFICIENT, true);
        aPoly[N - i].ModAddFastEq(one, Q);
        aPoly.SetFormat(Format::EVALUATION);
        m_monomials.push_back(std::move(aPoly));
    }
}
}
