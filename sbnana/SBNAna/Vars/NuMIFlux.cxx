#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana {

  //// ----------------------------------------------
  NuMIPpfxFluxWeight::NuMIPpfxFluxWeight()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    fFluxFilePath = std::string(sbndata) +
                   "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";

    TFile f(fFluxFilePath.c_str());
    if (f.IsZombie()) {
      std::cout << "NuMIPpfxFluxWeight: Failed to open " << fFluxFilePath << std::endl;
      std::abort();
    }

    for (int hcIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          std::string hNamePPFX = "ppfx_flux_weights/hweights_";
          if (hcIdx == 0)
            hNamePPFX += "fhc_";
          else
            hNamePPFX += "rhc_";
          if (flavIdx == 0)
            hNamePPFX += "nue";
          else
            hNamePPFX += "numu";
          if (signIdx == 1) hNamePPFX += "bar";

          TH1* h_ppfx = (TH1*)f.Get(hNamePPFX.c_str());
          if (!h_ppfx) {
            std::cout << "NuMIPpfxFluxWeight: failed to find " << hNamePPFX << " in " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }
  }

  NuMIPpfxFluxWeight::~NuMIPpfxFluxWeight()
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          delete fWeight[i][j][k];
  }

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeight([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxWeight(&slc->truth);
  });

  //// ----------------------------------------------

  const TruthVar kGetTruthNuMIFluxWeight([](const caf::SRTrueInteractionProxy* nu) -> double {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    if (!FluxWeightNuMI.fWeight[0][0][0]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;

    TH1* h = FluxWeightNuMI.fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(nu->E);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return 1.0;
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) return 1.0;
    return h->GetBinContent(bin);
  });

  /// --- VERSIONS WITH THE CONCRETE WEIGHTS
  //// ----------------------------------------------
  NuMIPpfxFluxWeightG3Chase::NuMIPpfxFluxWeightG3Chase()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIPpfxFluxG3ChaseWeight: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    fFluxFilePath = std::string(sbndata) +
                   "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";

    TFile f(fFluxFilePath.c_str());
    if (f.IsZombie()) {
      std::cout << "NuMIPpfxFluxG3ChaseWeight: Failed to open " << fFluxFilePath << std::endl;
      std::abort();
    }

    for (int hcIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          std::string hNamePPFX = "ppfx_flux_weights/hweights_";
          if (hcIdx == 0)
            hNamePPFX += "fhc_";
          else
            hNamePPFX += "rhc_";
          if (flavIdx == 0)
            hNamePPFX += "nue";
          else
            hNamePPFX += "numu";
          if (signIdx == 1) hNamePPFX += "bar";

          TH1* h_ppfx = (TH1*)f.Get(hNamePPFX.c_str());
          if (!h_ppfx) {
            std::cout << "NuMIPpfxFluxG3ChaseWeight: failed to find " << hNamePPFX << " in " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }

    fFluxFilePathG3Chase = std::string(sbndata) +
                           "beamData/NuMIdata/g3Chase_weights_rewritten.root";

    TFile f2(fFluxFilePathG3Chase.c_str());
    if (f2.IsZombie()) {
      std::cout << "NuMIPpfxFluxG3ChaseWeight: Failed to open " << fFluxFilePathG3Chase << std::endl;
      std::abort();
    }

    for (int flavIdx : {0, 1}) {
      for (int signIdx : {0, 1}) {
        for (int pdgIdx : {0, 1, 2, 3}) {

          if ( flavIdx==0 && pdgIdx==0 ) continue;

          std::string hNameG3Chase = "hweights_";
          hNameG3Chase += "fhc_";
          if (flavIdx == 0)
            hNameG3Chase += "nue";
          else
            hNameG3Chase += "numu";
          if (signIdx == 1) hNameG3Chase += "bar";
          if ( pdgIdx==0 ) hNameG3Chase += "_pipm";
          else if ( pdgIdx==1 ) hNameG3Chase += "_kpm";
          else if ( pdgIdx==2 ) hNameG3Chase += "_mu";
          else hNameG3Chase += "_k0l";

          TH1* h_g3chase = (TH1*)f2.Get(hNameG3Chase.c_str());
          if (!h_g3chase) {
            std::cout << "NuMIPpfxFluxG3ChaseWeight: failed to find " << hNameG3Chase << " in " << f2.GetName()
                      << std::endl;
            std::abort();
          }
          h_g3chase = (TH1*)h_g3chase->Clone(UniqueName().c_str());
          h_g3chase->SetDirectory(0);

          fWeightG3Chase[flavIdx][signIdx][pdgIdx] = h_g3chase;
        }
      }
    }
  }

  NuMIPpfxFluxWeightG3Chase::~NuMIPpfxFluxWeightG3Chase()
  {
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 4; ++k) {
          if ( k < 2 ) delete fWeight[i][j][k];
          delete fWeightG3Chase[i][j][k];
        }
      }
    }
  }

  unsigned int NuMIPpfxFluxWeightG3Chase::ParentPDGToIdx(int pdg) const
  {
    if      ( abs(pdg) == 211 ) return 0;
    else if ( abs(pdg) == 321 ) return 1;
    else if ( abs(pdg) == 13  ) return 2;
    else if ( abs(pdg) == 130 ) return 3;
    return 4;
  }


  //// ----------------------------------------------

  const Var kGetNuMIFluxWeightG3Chase([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxWeightG3Chase(&slc->truth);
  });

  //// ----------------------------------------------

  const TruthVar kGetTruthNuMIFluxWeightG3Chase([](const caf::SRTrueInteractionProxy* nu) -> double {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    /// Choose 1 1 1 for the G3Chase weight check since the 0 0 0 is meaningless here...
    if (!FluxWeightNuMIG3Chase.fWeight[0][0][0] || !FluxWeightNuMIG3Chase.fWeightG3Chase[1][1][1]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;
    unsigned int pdgIdx = FluxWeightNuMIG3Chase.ParentPDGToIdx(nu->parent_pdg);

    TH1* h = FluxWeightNuMIG3Chase.fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    double weight = 1.0;

    const int bin = h->FindBin(nu->E);
    if ( bin != 0 && bin != h->GetNbinsX() + 1 && !std::isinf(h->GetBinContent(bin)) && !std::isnan(h->GetBinContent(bin)) )
      weight*=h->GetBinContent(bin);

    if ( pdgIdx == 4 ) return weight;

    // return the weight as-is if looking for additional weight in the nue pion channel
    if ( pdgIdx == 0 && flavIdx == 0 ) return weight;

    TH1* h2 = FluxWeightNuMIG3Chase.fWeightG3Chase[flavIdx][signIdx][pdgIdx];
    assert(h2);

    const int bin2 = h2->FindBin(nu->E);
    if ( bin2 != 0 && bin2 != h2->GetNbinsX() + 1 && !std::isinf(h2->GetBinContent(bin2)) && !std::isnan(h2->GetBinContent(bin2)) )
      weight*=h2->GetBinContent(bin2);

    return weight;
  });

}
