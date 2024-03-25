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
  /// and as of 20 March 2024, updating to include an additional 
  /// weight for Kaon parents aimed at correcting for the G4 xsec error.
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

    // normal PPFX weights
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

    // Geometry updates, so-called "G3 chase"
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

    // Kaon correction for newer G4
    fFluxFilePathG4Kaon = std::string(sbndata) +
                          "beamData/NuMIdata/g4_10_4_weights_rewritten.root";

    TFile f3(fFluxFilePathG4Kaon.c_str());
    if (f3.IsZombie()) {
      std::cout << "NuMIPpfxFluxG3ChaseWeight: Failed to open " << fFluxFilePathG4Kaon << std::endl;
      std::abort();
    }

    for (int signIdx : {0, 1}) {
      std::string hNameKaon = "hweights_g4_10_4_kpm_fhc_numu";
      if (signIdx == 1) hNameKaon += "bar";

      TH1* h_g4kaon = (TH1*)f3.Get(hNameKaon.c_str());
      if (!h_g4kaon) {
        std::cout << "NuMIPpfxFluxG3ChaseWeight: failed to find " << hNameKaon << " in " << f3.GetName()
                  << std::endl;
        std::abort();
      }
      h_g4kaon = (TH1*)h_g4kaon->Clone(UniqueName().c_str());
      h_g4kaon->SetDirectory(0);

      fWeightG4Kaon[signIdx] = h_g4kaon;
    }

  }

  NuMIPpfxFluxWeightG3Chase::~NuMIPpfxFluxWeightG3Chase()
  {
    for (int i = 0; i < 2; ++i) {
      delete fWeightG4Kaon[i];
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 4; ++k) {
          if ( k < 2 ) delete fWeight[i][j][k];
          delete fWeightG3Chase[i][j][k];
        }
      }
    }
  }

  double NuMIPpfxFluxWeightG3Chase::GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu, const bool applyKaonRW) const
  {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    /// Choose 1 1 1 for the G3Chase weight check since the 0 0 0 is meaningless here...
    if (!fWeight[0][0][0] || !fWeightG3Chase[1][1][1]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;
    unsigned int pdgIdx = ParentPDGToIdx(nu->parent_pdg);

    TH1* h = fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    double weight = 1.0;

    const int bin = h->FindBin(nu->E);
    if ( bin != 0 && bin != h->GetNbinsX() + 1 && !std::isinf(h->GetBinContent(bin)) && !std::isnan(h->GetBinContent(bin)) )
      weight*=h->GetBinContent(bin);

    if ( pdgIdx == 4 ) return weight;

    // return the weight as-is if looking for additional weight in the nue pion channel
    if ( pdgIdx == 0 && flavIdx == 0 ) return weight;

    TH1* h2 = fWeightG3Chase[flavIdx][signIdx][pdgIdx];
    assert(h2);

    const int bin2 = h2->FindBin(nu->E);
    if ( bin2 != 0 && bin2 != h2->GetNbinsX() + 1 && !std::isinf(h2->GetBinContent(bin2)) && !std::isnan(h2->GetBinContent(bin2)) ) {
      weight*=h2->GetBinContent(bin2);
    }

    if ( pdgIdx==1 && applyKaonRW ) {
      TH1* h3 = fWeightG4Kaon[signIdx];
      assert(h3);
      const int bin3 = h3->FindBin(nu->E);
      if ( bin3 != 0 && bin3 != h3->GetNbinsX() + 1 && !std::isinf(h3->GetBinContent(bin3)) && !std::isnan(h3->GetBinContent(bin3)) ) {
        weight*=h3->GetBinContent(bin3);
      }
    }

    return weight;
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
    return FluxWeightNuMIG3Chase.GetWeightFromSRTrueInt(nu,false);
  });

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeightUpdated([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxWeightUpdated(&slc->truth);
  });

  //// ----------------------------------------------

  const TruthVar kGetTruthNuMIFluxWeightUpdated([](const caf::SRTrueInteractionProxy* nu) -> double {
    return FluxWeightNuMIG3Chase.GetWeightFromSRTrueInt(nu,true);
  });

}
