// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string>

class TH1;

namespace ana
{

  class NuMIPpfxFluxWeight
  {
  public:
    NuMIPpfxFluxWeight();
    ~NuMIPpfxFluxWeight();
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]

  protected:
    std::string fFluxFilePath;
  };

  // set up to use the flux weight
  static const NuMIPpfxFluxWeight FluxWeightNuMI;
  extern const Var kGetNuMIFluxWeight;
  extern const TruthVar kGetTruthNuMIFluxWeight;

  class NuMIPpfxFluxWeightG3Chase
  {
  public:
    // FOR NOW WE ONLY HAVE FHC...
    // parent idx --> pi+/-, K+/-, mu, K0L
    NuMIPpfxFluxWeightG3Chase();
    ~NuMIPpfxFluxWeightG3Chase();

    unsigned int ParentPDGToIdx(int pdg) const;

    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    mutable TH1* fWeightG3Chase[2][2][4]; // [nue/numu][nu/nubar][parent idx]

  protected:
    std::string fFluxFilePath;
    std::string fFluxFilePathG3Chase;
  };

  // set up to use the flux weight with G3Chase param controlling the concrete
  static const NuMIPpfxFluxWeightG3Chase FluxWeightNuMIG3Chase;
  extern const Var kGetNuMIFluxWeightG3Chase;
  extern const TruthVar kGetTruthNuMIFluxWeightG3Chase;

}
