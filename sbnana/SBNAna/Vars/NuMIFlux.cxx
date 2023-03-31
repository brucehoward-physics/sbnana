#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana
{

  //// ----------------------------------------------

  NuMIPpfxFluxWeight::NuMIPpfxFluxWeight()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if(!sbndata){
      std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup the sbndata product." << std::endl;
      std::abort();
    }

    const std::string fname = std::string(sbndata)+"/beamData/NuMIdata/icarus_numi_flux_syst_ana.root";
    TFile f(fname.c_str());
    if(f.IsZombie()){
      std::cout << "NuMIPpfxFluxWeight: Failed to open " << fname << std::endl;
      std::abort();
    }

    for(int hcIdx: {0, 1}){
      for(int flavIdx: {0, 1}){
        for(int signIdx: {0, 1}){
	  std::string      hNamePPFX = "flux_correction/htotal_";
          if(hcIdx == 0)   hNamePPFX += "fhc_"; else hNamePPFX += "rhc_";
          if(flavIdx == 0) hNamePPFX += "nue";  else hNamePPFX += "numu";
          if(signIdx == 1) hNamePPFX += "bar";

	  std::string      hNameNom = "ppfx_output/run15/";
          if(hcIdx == 0)   hNameNom += "fhc/nom/"; else hNameNom += "rhc/nom/";
          hNameNom+="hnom_";
          if(flavIdx == 0) hNameNom += "nue";      else hNameNom += "numu";
          if(signIdx == 1) hNameNom += "bar";

          TH1* h_ppfx = (TH1*)f.Get( hNamePPFX.c_str() );
          if(!h_ppfx){
	    std::cout << "NuMIPpfxFluxWeight: failed to find " << hNamePPFX << " in " << f.GetName() << std::endl;
	    std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          TH1* h_nom = (TH1*)f.Get( hNameNom.c_str() );
          if(!h_nom){
	    std::cout << "NuMIPpfxFluxWeight: failed to find " << hNameNom << " in " << f.GetName() << std::endl;
	    std::abort();
          }

          h_ppfx->Divide(h_nom);

          this->fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }
  }

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeight( [](const caf::SRSliceProxy* slc) -> double
	      {
		if ( slc->truth.index < 0 || abs(slc->truth.initpdg) == 16 ) return 1.0;

		if ( !FluxWeightNuMI.fWeight[0][0][0] ) {
		  std::cout << "Trying to access un-available weight array..." << std::endl;
		  std::abort();
		}

		unsigned int hcIdx = 0; // assume always FHC for now...
		unsigned int flavIdx = ( abs(slc->truth.initpdg) == 12 ) ? 0 : 1;
		unsigned int signIdx = ( slc->truth.initpdg > 0 ) ? 0 : 1;

		TH1* h = FluxWeightNuMI.fWeight[hcIdx][flavIdx][signIdx];
		assert(h);

		const int bin = h->FindBin( slc->truth.E );
		if(bin == 0 || bin == h->GetNbinsX()+1) return 1.0;
		return h->GetBinContent(bin);
	      });

  //// ----------------------------------------------

  NuMIPpfxFluxWeightNew::NuMIPpfxFluxWeightNew()
  {
    const std::string fname = "/pnfs/icarus/scratch/users/awood/numi_flux_systs_03042023/2023-03-04_out_450.37_7991.98_79512.66_qel7.root";
    TFile f(fname.c_str());
    if(f.IsZombie()){
      std::cout << "NuMIPpfxFluxWeightNew: Failed to open " << fname << std::endl;
      std::abort();
    }

    for(int hcIdx: {0, 1}){
      for(int flavIdx: {0, 1}){
        for(int signIdx: {0, 1}){
	  std::string      hName = "ppfx_flux_weights/hweights_";
	  if(hcIdx == 0)   hName += "fhc_"; else hName += "rhc_";
          if(flavIdx == 0) hName += "nue";  else hName += "numu";
          if(signIdx == 1) hName += "bar";

          TH1* h_ppfx = (TH1*)f.Get( hName.c_str() );
          if(!h_ppfx){
	    std::cout << "NuMIPpfxFluxWeightNew: failed to find " << hName << " in " << f.GetName() << std::endl;
	    std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          this->fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }
  }

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeightNew( [](const caf::SRSliceProxy* slc) -> double
	      {
		if ( slc->truth.index < 0 || abs(slc->truth.initpdg) == 16 ) return 1.0;

		if ( !FluxWeightNuMINew.fWeight[0][0][0] ) {
		  std::cout << "Trying to access un-available weight array..." << std::endl;
		  std::abort();
		}

		unsigned int hcIdx = 0; // assume always FHC for now...
		unsigned int flavIdx = ( abs(slc->truth.initpdg) == 12 ) ? 0 : 1;
		unsigned int signIdx = ( slc->truth.initpdg > 0 ) ? 0 : 1;

		TH1* h = FluxWeightNuMINew.fWeight[hcIdx][flavIdx][signIdx];
		assert(h);

		const int bin = h->FindBin( slc->truth.E );
		if(bin == 0 || bin == h->GetNbinsX()+1) return 1.0;
		return h->GetBinContent(bin);
	      });

  //// ----------------------------------------------

}
