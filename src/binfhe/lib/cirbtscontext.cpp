#include "cirbtscontext.h"
#include <cmath>
#include <unordered_map>
#include <cstdlib>
#include <iostream>

namespace lbcrypto{ 

void CirBTSContext::GenerateCirBTSContext(CirBTS_PARAMSET set, BINFHE_METHOD method) {
    constexpr double STD_DEV = 3.2;

    const std::unordered_map<CirBTS_PARAMSET, CirBTSContextParams> CircuitParamsMap({
            //                          numberBits|cyclOrder|latticeParam|  mod|   stdDev| BaseEP|  DigitsEP| BaseHT| DigitsHT| BaseSS| DigitsSS|BaseCC| DigitsCC|TraceShift|keyDist0|keyDist2 
        { STD128_CircuitBootstrap_CMUX_1, {54,      4096,      571,        1024,  STD_DEV,  1 << 26,   1,    1 << 17,   2,    1 << 28,    1,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        // Multi-bit PBS (k=2) experimental paramset; same security as CMUX_1, uses separate key path.
        { STD128_CircuitBootstrap_CMUX_1_MB2, {54,   4096,      571,        1024,  STD_DEV,  1 << 26,   1,    1 << 17,   2,    1 << 28,    1,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_1_EQ, {54,   4096,      571,        1024,  STD_DEV,  1 << 26,   1,    1 << 13,   3,    1 << 13,    3,    1 << 3,  4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_1_SUB1, {54, 4096,      571,        1024,  STD_DEV,  1 << 26,   1,    1 << 17,   2,    1 << 28,    1,    1 << 3,   4,   1,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_2, {54,      4096,      571,        1024,  STD_DEV,  1 << 17,   2,    1 << 17,   2,    1 << 28,    1,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_3, {54,      4096,      571,        1024,  STD_DEV,  1 << 17,   2,    1 << 17,   2,    1 << 19,    2,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_4, {54,      4096,      571,        1024,  STD_DEV,  1 << 17,   2,    1 << 13,   3,    1 << 19,    2,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
        { STD128_CircuitBootstrap_CMUX_5, {54,      4096,      571,        1024,  STD_DEV,  1 << 17,   2,    1 << 11,   4,    1 << 19,    2,    1 << 3,   4,   0,   UNIFORM_BINARY, UNIFORM_BINARY} },
    });

    auto search = CircuitParamsMap.find(set);
    if (CircuitParamsMap.end() == search ){
        std::string errMsg("ERROR: Unknown parameter set [" + std::to_string(set) + "] for circuitbootstrap");
        OPENFHE_THROW(config_error, errMsg);
    }

    CirBTSContextParams params = search->second;
    const bool verbose_overrides = std::getenv("CIRBTS_PARAM_VERBOSE") != nullptr;
    auto override_usint = [&](const char* env_name, usint& target) {
        const char* v = std::getenv(env_name);
        if (!v || !*v) {
            return;
        }
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(v, &end, 10);
        if (end != v && end && *end == '\0') {
            target = static_cast<usint>(parsed);
            if (verbose_overrides) {
                std::cout << "[CIRBTS] Override " << env_name << "=" << target << std::endl;
            }
        }
    };
    override_usint("CIRBTS_BASE_EP", params.BaseEP);
    override_usint("CIRBTS_DIGITS_EP", params.DigitsEP);
    override_usint("CIRBTS_BASE_HT", params.BaseHT);
    override_usint("CIRBTS_DIGITS_HT", params.DigitsHT);
    override_usint("CIRBTS_BASE_SS", params.BaseSS);
    override_usint("CIRBTS_DIGITS_SS", params.DigitsSS);
    override_usint("CIRBTS_BASE_CC", params.BaseCC);
    override_usint("CIRBTS_DIGITS_CC", params.DigitsCC);
    override_usint("CIRBTS_TRACE_SHIFT", params.TraceShift);

    //level 2 prime modulus 
    NativeInteger Q(LastPrime<NativeInteger>(params.numberBits, params.cyclOrder));

    usint ringDim = params.cyclOrder / 2;
    auto lweparams = std::make_shared<LWECryptoParams>(params.latticeParam, ringDim, params.mod, Q, params.mod,
                                                        params.stdDev, 1, params.keyDist0);
    uint32_t baseR_ep = params.BaseEP;
    uint32_t baseR_cc = params.BaseCC;
    auto rgswparams1 = std::make_shared<RingGSWCryptoParams>(ringDim, Q, params.mod, params.BaseEP, baseR_ep,
                                                             method, params.stdDev, params.DigitsEP, params.keyDist2, false, 10);
    auto rlweparams = std::make_shared<RLWECryptoParams>(params.cyclOrder, ringDim, Q, params.stdDev, params.BaseHT,
                                                        params.DigitsHT, params.BaseSS, params.DigitsSS, params.keyDist2);
    if (params.TraceShift != 0) {
        const uint32_t logN = static_cast<uint32_t>(std::log2(static_cast<double>(ringDim)));
        if (params.TraceShift >= logN) {
            OPENFHE_THROW(config_error, "TraceShift must be smaller than log2(N)");
        }
        rlweparams->SetTraceShift(params.TraceShift);
    }
    auto rgswparams2 = std::make_shared<RingGSWCryptoParams>(ringDim, Q, params.mod, params.BaseCC, baseR_cc,
                                                             method, params.stdDev, params.DigitsCC, params.keyDist2, false, 10);
    m_params = std::make_shared<CirBTSCryptoParams>(lweparams, rgswparams1, rlweparams, rgswparams2);
    m_cirbtsscheme = std::make_shared<CirBTSScheme>(method);
}

RLWEPrivateKey CirBTSContext::RLWEKeyGen() const{
    auto& RLWEParams = m_params->GetRLWEParams();
    if (RLWEParams->GetKeyDist() == GAUSSIAN)
        return m_RLWEscheme->KeyGenGaussian(RLWEParams->GetN(), RLWEParams->GetQ(), 3.2);
    else if (RLWEParams->GetKeyDist() == UNIFORM_TERNARY)
    {
        return m_RLWEscheme->KeyGenTernary(RLWEParams->GetN(), RLWEParams->GetQ());
    }
    return m_RLWEscheme->KeyGenBinary(RLWEParams->GetN(), RLWEParams->GetQ());
}

LWEPrivateKey CirBTSContext::KeyGen() const{
    auto& LWEParams = m_params->GetLWEParams();
    if (LWEParams->GetKeyDist() == GAUSSIAN)
        return m_LWEscheme->KeyGenGaussian(LWEParams->Getn(), LWEParams->GetqKS());
    else if (LWEParams->GetKeyDist() == UNIFORM_TERNARY)
    {
        return m_LWEscheme->KeyGen(LWEParams->Getn(), LWEParams->GetqKS());
    }
    return m_LWEscheme->KeyGenBinary(LWEParams->Getn(), LWEParams->GetqKS());
}

LWECiphertext CirBTSContext::Encrypt(ConstLWEPrivateKey& sk, LWEPlaintext m, LWEPlaintextModulus p) const{
    const auto& LWEParams = m_params->GetLWEParams();
    LWECiphertext ct = m_LWEscheme->Encrypt(LWEParams, sk, m, p, LWEParams->Getq());
    return ct;
}

void CirBTSContext::Decrypt(ConstLWEPrivateKey& sk, ConstLWECiphertext& ct, LWEPlaintext* result, LWEPlaintextModulus p) const{
    auto&& LWEParams = m_params->GetLWEParams();
    m_LWEscheme->Decrypt(LWEParams, sk, ct, result, p);
}

void CirBTSContext::CirBTKeyGen(ConstLWEPrivateKey& sk, ConstRLWEPrivateKey skNTT, KEYGEN_MODE keygenMode){
    m_BTKey           = m_cirbtsscheme->KeyGen(m_params, sk, skNTT, keygenMode);
}

RGSWCiphertext CirBTSContext::CircuitBootstrapping(ConstLWECiphertext& ct) const{
    return m_cirbtsscheme->CircuitBootstrap(m_params, m_BTKey, ct);
}

}
