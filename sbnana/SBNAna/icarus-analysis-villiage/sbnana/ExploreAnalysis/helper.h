// Place to put our Cuts and Vars

// SBNAna includes
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

//#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

//#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
//#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

#include "TGraphAsymmErrors.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TPaveText.h"

// C++
#include <vector>
#include <string>

using namespace ana;

/////////////

// First attempts at NuMI nu direction vectors
TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);
// also
TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

/////////////
// Plot stylings
/////////////

void PrintICARUSSim()
{
  TPaveText *tpt = new TPaveText(0.63,0.91,0.9,0.96,"NB NDC");
  TText *txt = tpt->AddText("ICARUS Simulation");
  txt->SetTextColor(kGray+1);

  tpt->Draw();
}

void PrintICARUSData()
{
  TPaveText *tpt = new TPaveText(0.63,0.91,0.9,0.96,"NB NDC");
  TText *txt = tpt->AddText("ICARUS Data");
  txt->SetTextColor(kRed+1);

  tpt->Draw();
}

void PrintPreliminary()
{
  TPaveText *tpt = new TPaveText(0.25,0.91,0.52,0.96,"NB NDC");
  TText *txt = tpt->AddText("Preliminary");
  txt->SetTextColor(kBlue+1);

  tpt->Draw();
}

/////////////////////////////////////////////////
// Function to make an error band as in
// EnsembleSpectrum, but taking 3 spectra:
//   Nominal, +1 sigma, -1 sigma
// See https://github.com/SBNSoftware/sbnana/blob/develop/sbnana/CAFAna/Core/EnsembleSpectrum.cxx
/////////////////////////////////////////////////
TGraphAsymmErrors* ErrorBandFromSpectra( const Spectrum &sNom, const Spectrum &sMinus, const Spectrum &sPlus,
					 double exposure, EExposureType expotype=kPOT, EBinType bintype=kBinContent)
{
  std::unique_ptr<TH1D> hnom(sNom.ToTH1(exposure, expotype, bintype));
  std::unique_ptr<TH1D> hplus(sPlus.ToTH1(exposure, expotype, bintype));
  std::unique_ptr<TH1D> hminus(sMinus.ToTH1(exposure, expotype, bintype));

  TGraphAsymmErrors* g = new TGraphAsymmErrors;

  for(int binIdx = 0; binIdx < hnom->GetNbinsX()+2; ++binIdx){
    const double xnom = hnom->GetXaxis()->GetBinCenter(binIdx);
    const double ynom = hnom->GetBinContent(binIdx);
    g->SetPoint(binIdx, xnom, ynom);

    const double dx = hnom->GetXaxis()->GetBinWidth(binIdx);

    // BH: I think for *this* particular functionality, I actually want to take the max and min of the + and - ...
    const double y0 = std::min(hminus->GetBinContent(binIdx), hplus->GetBinContent(binIdx));
    const double y1 = std::max(hplus->GetBinContent(binIdx), hminus->GetBinContent(binIdx));

    // It's theoretically possible for the central value to be outside the
    // error bands - clamp to zero in that case
    g->SetPointError(binIdx, dx/2, dx/2,
		     std::max(ynom-y0, 0.),
		     std::max(y1-ynom, 0.));
  }

  return g;
}

// Also let's make a helper function for summing in quadrature... taking also the nominal to set the x and y
TGraphAsymmErrors* ErrorBandInQuad( const Spectrum &sNom, const std::vector<TGraphAsymmErrors*> &errGraphs,
				    double exposure, EExposureType expotype=kPOT, EBinType bintype=kBinContent)
{
  std::unique_ptr<TH1D> hnom(sNom.ToTH1(exposure, expotype, bintype));

  TGraphAsymmErrors* g = new TGraphAsymmErrors;

  for(int binIdx = 0; binIdx < hnom->GetNbinsX()+2; ++binIdx){
    const double xnom = hnom->GetXaxis()->GetBinCenter(binIdx);
    const double ynom = hnom->GetBinContent(binIdx);
    g->SetPoint(binIdx, xnom, ynom);

    const double dx = hnom->GetXaxis()->GetBinWidth(binIdx);

    double sumHi=0.;
    double sumLo=0.;

    for(auto const& hErr : errGraphs){
      sumHi+=std::pow(hErr->GetErrorYhigh(binIdx),2);
      sumLo+=std::pow(hErr->GetErrorYlow(binIdx),2);
    }

    double errHi = std::sqrt(sumHi);
    double errLo = std::sqrt(sumLo);

    std::cout << binIdx << " high = " << errHi << ", low = " << errLo << std::endl;

    g->SetPointError(binIdx, dx/2, dx/2,
		     errHi,
		     errLo);
  }

  return g;
}

/////////////////////////////////////////////////
// Colors
/////////////////////////////////////////////////

// We'll use it as Signal, Other NuMu CC, NC, Cosmics ('out-of-time'), In-time cosmic
int colorwheel[5] = {kBlue, kGreen+2, kRed, kMagenta, kViolet-1};

// order: QE, MEC, RES, DIS, COH, NonSigCC
int colorwheel_mode[6] = {kBlue, kCyan+1, kGreen+2, kOrange+1, kGreen, kGray+2};

// order: true signal, numu mu wrong, numu p wrong, numu both wrong, numu non fid, numu other, nu other, cosmic, cosmic in-time
int colorwheel_class[9] = {kBlue, kAzure+7, kRed, kGray+2, kGreen+2, kGreen, kOrange+2, kMagenta, kViolet-1 };

/////////////////////////////////////////////////
// Binning
/////////////////////////////////////////////////

const Binning kBinsE = Binning::Simple(50,0.,5.);
//const Binning kBinsP = Binning::Simple(25,0.,2.5);
const Binning kBinsP = Binning::Simple(35,0.,3.5);

const Binning kBinsPZoom = Binning::Simple(25,0.,0.5);

const Binning kBinsProtonPZoom = Binning::Simple(50,0.,1.5);

const Binning kBinsProtonP = Binning::Simple(16,0,1.6);
const Binning kBinsCounts = Binning::Simple(8,-0.5,7.5);

const Binning kBinsT = Binning::Simple(30,0,1.5); // |t| bins
const Binning kBinsTzoom = Binning::Simple(20,0,.2);

const Binning kBinsPID = Binning::Simple(30,0,150);
const Binning kBinsPIDA = Binning::Simple(25,0,50);

const Binning kBinsL = Binning::Simple(50,0.,1000.);
const Binning kBinsLE = Binning::Simple(40,0.,4.);

const Binning kBinsCosTh = Binning::Simple(10,-1.,1.);

const Binning kBinsRTLen = Binning::Simple(31,-0.05,3.05);

const Binning kBinsStubLen = Binning::Simple(10,0.,5.);
const Binning kBinsStubDQDx = Binning::Simple(16,0.,8.e5);

const Binning kBinsPDGType = Binning::Simple(8,0,8);

const Binning kBinsP_Custom = Binning::Custom({0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.5, 3.5});

/////////////////////////////////////////////////\/////////////////////////////////////////////////
/////////////////////////////////////////////////\/////////////////////////////////////////////////

// SLICE LEVEL
bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
	    ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
	  (( y > -181.86 + 25 && y < 134.96 - 25 ) &&
	   ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

const Cut isInFV_Vtx([](const caf::SRSliceProxy* sr)
		     {
		       const auto& vtx = sr->truth.position;

		       if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

		       return (( ( vtx.x < -61.94 - 25 && vtx.x > -358.49 + 25 ) ||
				 ( vtx.x >  61.94 + 25 && vtx.x <  358.49 - 25 )) &&
			       (( vtx.y > -181.86 + 25 && vtx.y < 134.96 - 25 ) &&
				( vtx.z > -894.95 + 30 && vtx.z < 894.95 - 50 ) ));
		     });

const Cut isInAV_Vtx([](const caf::SRSliceProxy* sr)
                     {
                       const auto& vtx = sr->truth.position;

                       if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

                       return (( ( vtx.x < -61.94 && vtx.x > -358.49 ) ||
                                 ( vtx.x >  61.94 && vtx.x <  358.49 )) &&
                               (( vtx.y > -181.86 && vtx.y < 134.96 ) &&
                                ( vtx.z > -894.95 && vtx.z < 894.95 ) ));
                     });

/////////////////////////////////////////////////
// Signal definition(s)
/////////////////////////////////////////////////

/*
  -- USE THESE FROM NumuCutsIcarus202106 if needed --
  const Cut kIsNuSlice = ( kTruthIndex >= 0.f );
  const Cut kIsCosmic = ( !kIsNuSlice );
*/

const Cut kCutTrueSig([](const caf::SRSliceProxy* slc)
		      {
			if(slc->truth.index < 0) return false;

			return ( abs(slc->truth.pdg) == 14 &&
				 slc->truth.iscc &&
				 !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
				 isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) );
		      });

const Cut kCutTrueSig_50cmMuCont_100cmMuExit([](const caf::SRSliceProxy* slc)
	{
  	if(slc->truth.index < 0) return false;

    bool signalLepton = false;

    if ( abs(slc->truth.pdg) == 14 &&
				 slc->truth.iscc &&
				 !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
				 isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) ) {
      for ( auto const& prim : slc->truth.prim ) {
        if ( prim.pdg == 2212 && prim.length > 100. ) {
          signalLepton = true;
          break;
        }
        else if ( prim.pdg == 2212 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          break;
        }
      }
    }

		return signalLepton;
  });

const Cut kCutTrueSigOnAr([](const caf::SRSliceProxy* slc)
			  {
			    if(slc->truth.index < 0) return false;

			    return ( abs(slc->truth.pdg) == 14 &&
				     slc->truth.targetPDG == 1000180400 &&
				     slc->truth.iscc &&
				     !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
				     isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) );
			  });

const Cut kCutTrueSigContained([](const caf::SRSliceProxy* slc)
                      {
			if(slc->truth.index < 0) return false;

			if( abs(slc->truth.pdg) == 14 &&
			    slc->truth.iscc &&
			    !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
			    isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) )
			  {
			    for ( auto const& iprim : slc->truth.prim ) {
			      if ( iprim.pdg == 13 && iprim.contained ) return true;
			    }
			  }

			return false;
                      });

const Cut kCutNuCCButNotSigAll([](const caf::SRSliceProxy* slc)
                      {
                        if(slc->truth.index < 0) return false;

                        bool isSignal = (abs(slc->truth.pdg) == 14 &&
                            slc->truth.iscc &&
                            !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                            isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));

                        return ( slc->truth.iscc && !isSignal );
                      });

const Cut kCutNuCCButNotSigContained([](const caf::SRSliceProxy* slc)
				     {
				       if(slc->truth.index < 0) return false;

				       bool isSignal = false;
				       if( abs(slc->truth.pdg) == 14 &&
					   slc->truth.iscc &&
					   !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
					   isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) )
					 {
					   for ( auto const& iprim : slc->truth.prim ) {
					     if ( iprim.pdg == 13 && iprim.contained ) {
					       isSignal = true;
					       break;
					     }
					   }
					 }

				       return ( slc->truth.iscc && !isSignal );
				     });

const Cut kCutTrueNC([](const caf::SRSliceProxy* slc)
		 {
		   if(slc->truth.index < 0) return false;
		   return !slc->truth.iscc;
		 });

const Cut kCutCosmic([](const caf::SRSliceProxy* slc)
		     {
		       return (slc->truth.index < 0);
		     });

// kCutTrueSig && Interaction Type
const Cut kCutTrueSigQEL([](const caf::SRSliceProxy* slc)
			 {
			   if( !kCutTrueSig(slc) ) return false;

			   return slc->truth.genie_mode==caf::kQE && slc->truth.genie_inttype==caf::kCCQE;
			 });

const Cut kCutTrueSigMEC([](const caf::SRSliceProxy* slc)
                         {
                           if( !kCutTrueSig(slc) ) return false;

                           return slc->truth.genie_mode==caf::kMEC;
                         });

const Cut kCutTrueSigRES([](const caf::SRSliceProxy* slc)
                         {
                           if( !kCutTrueSig(slc) ) return false;

                           return slc->truth.genie_mode==caf::kRes;
                         });

const Cut kCutTrueSigDIS([](const caf::SRSliceProxy* slc)
                         {
                           if( !kCutTrueSig(slc) ) return false;

                           return slc->truth.genie_mode==caf::kDIS;
                         });

const Cut kCutTrueSigCOH([](const caf::SRSliceProxy* slc)
                         {
                           if( !kCutTrueSig(slc) ) return false;

                           return slc->truth.genie_mode==caf::kCoh;
                         });

const Cut kCutTrueSigELS([](const caf::SRSliceProxy* slc)
                         {
                           if( !kCutTrueSig(slc) ) return false;

                           return !( (slc->truth.genie_mode==caf::kQE && slc->truth.genie_inttype==caf::kCCQE) ||
                                      slc->truth.genie_mode==caf::kMEC ||
                                      slc->truth.genie_mode==caf::kRes ||
                                      slc->truth.genie_mode==caf::kDIS ||
                                      slc->truth.genie_mode==caf::kCoh );
                         });


// 1 mu + >=1 p with 0 charged pi
const Cut kCutSignal_1muNp0pi([](const caf::SRSliceProxy* slc)
		      {
      			if(slc->truth.index < 0) return false;

      			bool signal = (abs(slc->truth.pdg) == 14 &&
                           slc->truth.iscc &&
                           !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                           isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));

            if ( !signal ) return false;

            unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
            for ( auto const& prim : slc->truth.prim ) {
              if ( abs(prim.pdg) == 13 ) nMu+=1;
              if ( abs(prim.pdg) == 2212 ) nP+=1;
              if ( abs(prim.pdg) == 111 ) nPi0+=1;
              if ( abs(prim.pdg) == 211 ) nChgPi+=1;
            }

            return (nMu==1 && nP>=1 && nPi0 == 0 && nChgPi == 0);
		      });

const Cut kCutNuCC_Not_1muNp0pi([](const caf::SRSliceProxy* slc)
		      {
      			if(slc->truth.index < 0) return false;
            if(!slc->truth.iscc ) return false;

      			bool signal = (abs(slc->truth.pdg) == 14 &&
                           slc->truth.iscc &&
                           !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                           isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));

            unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
            for ( auto const& prim : slc->truth.prim ) {
              if ( abs(prim.pdg) == 13 ) nMu+=1;
              if ( abs(prim.pdg) == 2212 ) nP+=1;
              if ( abs(prim.pdg) == 111 ) nPi0+=1;
              if ( abs(prim.pdg) == 211 ) nChgPi+=1;
            }

            bool topology = (nMu==1 && nP>=1 && nPi0 == 0 && nChgPi == 0);

            if ( !signal || !topology ) return true;
            return false;
		      });



// --> What non signal selections should I look at?
// ----> NumuCC Muon is < 50cm (contained) or 100cm (uncontained) but reco as >
// ----> NumuCC Proton is contained (reco & true) but P < 400 MeV/c in truth, or Proton is uncontained but reco as contained --> grouping together as "w/o proton"
// ----> NumuCC non-fiducial
// ----> NumuCC other
// ----> Nu other? (e.g. NC with a charged pi)
// ----> Cosmics

// Do these with slice cuts on the slice truth (as above)
const Cut k1Mu1P_TrueSigTopology([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 &&
                                             slc->truth.iscc &&
                                             !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                                             isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));
                            if ( !isSignal ) return false;

                            // Check muon length and proton momentum
                            bool signalLepton = false;
                            bool signalHadron = false;

                            for ( auto const& prim : slc->truth.prim ) {
                              if ( signalLepton && signalHadron ) break;
                              if ( prim.pdg == 13 && prim.length > 100. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 13 && prim.length > 50. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 2212 && !signalHadron ) {
                                if ( !prim.contained ) continue;
                                float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
                                if ( p > 0.4 )
                                  signalHadron = true;
                                continue;
                              }
                            }

                            if ( !signalLepton || !signalHadron ) return false;
                            return true;
                          });

const Cut k1Mu1P_SliceIsNumuCC_MuLengthWrong([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 &&
                                             slc->truth.iscc &&
                                             !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                                             isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));
                            if ( !isSignal ) return false;

                            // Check muon length and proton momentum
                            bool signalLepton = false;
                            bool signalHadron = false;

                            for ( auto const& prim : slc->truth.prim ) {
                              if ( signalLepton && signalHadron ) break;
                              if ( prim.pdg == 13 && prim.length > 100. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 13 && prim.length > 50. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 2212 && !signalHadron ) {
                                if ( !prim.contained ) continue;
                                float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
                                if ( p > 0.4 )
                                  signalHadron = true;
                                continue;
                              }
                            }

                            if ( !signalLepton && signalHadron ) return true;
                            return false;
                          });

const Cut k1Mu1P_SliceIsNumuCC_PWrong([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 &&
                                             slc->truth.iscc &&
                                             !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                                             isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));
                            if ( !isSignal ) return false;

                            // Check muon length and proton momentum
                            bool signalLepton = false;
                            bool signalHadron = false;

                            for ( auto const& prim : slc->truth.prim ) {
                              if ( signalLepton && signalHadron ) break;
                              if ( prim.pdg == 13 && prim.length > 100. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 13 && prim.length > 50. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 2212 && !signalHadron ) {
                                if ( !prim.contained ) continue;
                                float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
                                if ( p > 0.4 )
                                  signalHadron = true;
                                continue;
                              }
                            }

                            if ( signalLepton && !signalHadron ) return true;
                            return false;
                          });

const Cut k1Mu1P_SliceIsNumuCC_BothWrong([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 &&
                                             slc->truth.iscc &&
                                             !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                                             isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));
                            if ( !isSignal ) return false;

                            // Check muon length and proton momentum
                            bool signalLepton = false;
                            bool signalHadron = false;

                            for ( auto const& prim : slc->truth.prim ) {
                              if ( signalLepton && signalHadron ) break;
                              if ( prim.pdg == 13 && prim.length > 100. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 13 && prim.length > 50. && !signalLepton ) {
                                signalLepton = true;
                                continue;
                              }
                              else if ( prim.pdg == 2212 && !signalHadron ) {
                                if ( !prim.contained ) continue;
                                float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
                                if ( p > 0.4 )
                                  signalHadron = true;
                                continue;
                              }
                            }

                            if ( !signalLepton && !signalHadron ) return true;
                            return false;
                          });

const Cut k1Mu1P_SliceIsNumuCC_NonFiducial([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 &&
                                             slc->truth.iscc &&
                                             !std::isnan(slc->truth.position.x) && !std::isnan(slc->truth.position.y) && !std::isnan(slc->truth.position.z) &&
                                             !isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z));
                            if ( !isSignal ) return false;

                            return true;
                          });

const Cut k1Mu1P_SliceIsNumuCC_GENERAL_FOR_OTHER([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            bool isSignal = (abs(slc->truth.pdg) == 14 && slc->truth.iscc);
                            if ( !isSignal ) return false;

                            return true;
                          });
const Cut k1Mu1P_SliceIsNumuCC_Other = k1Mu1P_SliceIsNumuCC_GENERAL_FOR_OTHER && 
                                       !k1Mu1P_SliceIsNumuCC_NonFiducial &&
                                       !k1Mu1P_SliceIsNumuCC_BothWrong &&
                                       !k1Mu1P_SliceIsNumuCC_PWrong &&
                                       !k1Mu1P_SliceIsNumuCC_MuLengthWrong &&
                                       !k1Mu1P_TrueSigTopology;

const Cut k1Mu1P_SliceIsNuOther([](const caf::SRSliceProxy* slc)
                          {
                            if(slc->truth.index < 0) return false;

                            if (abs(slc->truth.pdg) == 14 && slc->truth.iscc) return false;
                            return true;
                          });

const Cut k1Mu1P_SliceIsCosmic([](const caf::SRSliceProxy* slc)
                          {
                            return slc->truth.index < 0;
                          });


// kCutTrueSig && Final State?


// To pick "good" matches
/*
const Cut kSliceEffCut([](const caf::SRSliceProxy* slc)
		       {
			 if(slc->truth.index < 0) return false;

			 return slc->tmatch.eff > 0.5;
		       });
*/

/////////////////////////////////////////////////
// Slice Vars
/////////////////////////////////////////////////

const Var kTrueMuonMomentum([](const caf::SRSliceProxy* slc) -> float
			    {
			      if ( abs(slc->truth.pdg) != 14 ) return -5.f;

			      bool foundPdg=false;
			      float startP(0.);
			      for ( auto const& iPrim : slc->truth.prim ) {
				if ( abs(iPrim.pdg)==13 ) {
				  startP = std::sqrt( std::pow(iPrim.genp.x,2) +
						      std::pow(iPrim.genp.y,2) +
						      std::pow(iPrim.genp.z,2) );
				  foundPdg=true;
				  break;
				}
			      }
			      if ( !foundPdg ) return -5.f;

			      return startP;
			    });

/////////////////////////////////////////////////
// Slice Vars and Cuts (that aren't already in another header or in the fiducial volume stuff above)
/////////////////////////////////////////////////


// Redo the following two slice vars from NumuVars to use 10cm for containment (and to consider the other cryostat)

const Var kPTrackIndNew([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.trk.size(); ++i)
    {
        auto const& trk = slc->reco.trk.at(i);
        if(trk.bestplane == -1) continue;

        // First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && trk.pfp.parent_is_primary);

        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

        const bool Contained = (!isnan(trk.end.x) &&
                                ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                 (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                !isnan(trk.end.y) &&
                                ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                !isnan(trk.end.z) &&
                                ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
    	  {
	        Longest = trk.len;
	        PTrackInd = i;
    	  }
    }
    return PTrackInd;
  });

const Var kRecoMuonBestMatchPDG([](const caf::SRSliceProxy* slc) -> float {
    float pdg(-1.f);

    if ( kPTrackIndNew(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNew(slc));
      if ( abs(trk.truth.p.pdg) == 11 )        pdg = 0.;
      else if ( abs(trk.truth.p.pdg) == 13 )   pdg = 1.;
      else if ( abs(trk.truth.p.pdg) == 22 )   pdg = 2.;
      else if ( abs(trk.truth.p.pdg) == 111 )  pdg = 3.;
      else if ( trk.truth.p.pdg == 211 )       pdg = 4.;
      else if ( trk.truth.p.pdg == -211 )      pdg = 5.;
      else if ( abs(trk.truth.p.pdg) == 2212 ) pdg = 6.;
      else                                     pdg = 7.;
    }

    return pdg + std::numeric_limits<float>::epsilon();
  });

const MultiVar kTrueParticlePionPDGs([](const caf::SRSliceProxy* slc) {
    std::vector<double> pdgs;

    if ( slc->truth.index < 0 ) return pdgs;

    for ( auto const& prim : slc->truth.prim ) {
      if ( abs(prim.pdg) == 22 ) pdgs.push_back( 2. + std::numeric_limits<double>::epsilon() );
      if ( abs(prim.pdg) == 111 ) pdgs.push_back( 3. + std::numeric_limits<double>::epsilon() );
      if ( prim.pdg == 211 ) pdgs.push_back( 4. + std::numeric_limits<double>::epsilon() );
      if ( prim.pdg == -211 ) pdgs.push_back( 5. + std::numeric_limits<double>::epsilon() );
    }

    if ( pdgs.size() == 0 ) pdgs.push_back( 7. + std::numeric_limits<double>::epsilon() );

    return pdgs;
  });

const Var kRecoMuonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNew(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackIndNew(slc));
        const bool Contained = ( !isnan(trk.end.x) &&
				                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				                          (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) p = trk.rangeP.p_muon;
        else p = trk.mcsP.fwdP_muon;
      }
    return p;
  });

const Cut kRecoMuonContained([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNew(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackIndNew(slc));
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                  (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) return true;
        else return false;
      }
    return false;
  });

// PROTONs
// -- just picking the best proton track
const Var kPTrackIndNewProton([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kPTrackIndNew(slc);

    int idxScdy = -1;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; // need a different track...
  
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      float angle = -5.0;
      if ( primaryInd >= 0 ) {
        const unsigned int idxPrim = (unsigned int)primaryInd;
        TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
        TVector3 pDir( slc->reco.trk[idxTrk].dir.x, slc->reco.trk[idxTrk].dir.y, slc->reco.trk[idxTrk].dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      // do we want to make the proton cut even tighter on PID
      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = (int)idxTrk;
      }
    }

    return idxScdy;
  }); // kPTrackIndNewProton

const Var kPTrackIndNewProtonSimple([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kPTrackIndNew(slc);

    int idxScdy = -1;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; // need a different track...
  
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      float angle = -5.0;
      if ( primaryInd >= 0 ) {
	      const unsigned int idxPrim = (unsigned int)primaryInd;
	      TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
	      TVector3 pDir( slc->reco.trk[idxTrk].dir.x, slc->reco.trk[idxTrk].dir.y, slc->reco.trk[idxTrk].dir.z );
	      angle = TMath::Cos(muDir.Angle(pDir));
      }

      // do we want to make the proton cut even tighter on PID
      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && angle >= -0.9 && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = (int)idxTrk;
      }
    }

    return idxScdy;
  }); // kPTrackIndNewProtonSimple

const Var kRecoProtonP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNewProton(slc) >= 0 )
    {
        auto const& trk = slc->reco.trk.at(kPTrackIndNewProton(slc));
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                  (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) p = trk.rangeP.p_proton;
        else {
	        std::cout << "Currently kPTrackIndNewProton requires a contained proton... Why am I trying to use MCS here??" << std::endl;
	        p = trk.mcsP.fwdP_proton;
      	}
    }
    return p;
  });

const Var kRecoProtonPSimple([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                  (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) p = trk.rangeP.p_proton;
        else {
	        std::cout << "Currently kPTrackIndNewProton requires a contained proton... Why am I trying to use MCS here??" << std::endl;
          p = trk.mcsP.fwdP_proton;
        }
      }
    return p;
  });

const Var kRecoProtonChi2PSimple([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));
      const bool Contained = ( !isnan(trk.end.x) &&
			       ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			       !isnan(trk.end.y) &&
			       ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			       !isnan(trk.end.z) &&
			       ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2proton = trk.chi2pid[trk.bestplane].chi2_proton;
    }

    return chi2proton;
  });

const Var kRecoProtonChi2MuSimple([](const caf::SRSliceProxy* slc) -> float {
    float chi2muon(-5.f);

    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));
      const bool Contained = ( !isnan(trk.end.x) &&
			       ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			       !isnan(trk.end.y) &&
			       ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			       !isnan(trk.end.z) &&
			       ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2muon = trk.chi2pid[trk.bestplane].chi2_muon;
    }

    return chi2muon;
  });

// PIDA

const Var kRecoProtonPIDASimple([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));
      const bool Contained = ( !isnan(trk.end.x) &&
			       ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			       !isnan(trk.end.y) &&
			       ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			       !isnan(trk.end.z) &&
			       ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2proton = trk.chi2pid[trk.bestplane].pida;
    }

    return chi2proton;
  });

const Var kRecoMuonPIDA([](const caf::SRSliceProxy* slc) -> float {
    float chi2muon(-5.f);

    if ( kPTrackIndNew(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNew(slc));
      const bool Contained = ( !isnan(trk.end.x) &&
			       ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			       !isnan(trk.end.y) &&
			       ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			       !isnan(trk.end.z) &&
			       ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2muon = trk.chi2pid[trk.bestplane].pida;
    }

    return chi2muon;
  });


// Is selected proton (kPTrackIndNewProton) actually a proton
const Cut kRecoProtonIsTrueProton([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProton(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProton(slc));

      if ( trk.truth.p.pdg == 2212 ) return true;
      else return false;
    }

    return false;
  });

const Cut kRecoProtonIsTrueProtonSimple([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));

      if ( trk.truth.p.pdg == 2212 ) return true;
      else return false;
    }

    return false;
  });

const Cut kRecoProtonIsTrueStoppingProtonSimple([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProtonSimple(slc));

      if ( trk.truth.p.pdg == 2212 && trk.truth.p.daughters.size() == 0 ) return true;
      else return false;
    }

    return false;
  });

const Var kMatchedProtonTrueMomentum([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNewProton(slc) >= 0 )
    {
      auto const& trk = slc->reco.trk.at(kPTrackIndNewProton(slc));

      if ( trk.truth.p.pdg == 2212 ) p = sqrt( std::pow( trk.truth.p.startp.x, 2 ) + std::pow( trk.truth.p.startp.y, 2 ) + std::pow( trk.truth.p.startp.z, 2 ) );
    }

    return p;
  });

// -- multiplicities

// THIS VERSION DOESN'T SEEM TO WORK?
const Var kProtonMult_True_PreFSI([](const caf::SRSliceProxy* slc) -> float {
    if ( slc->truth.index < 0 ) return -5.0f;

    return (float)slc->truth.nproton;
  });

const Var kProtonMult_True_PostFSI([](const caf::SRSliceProxy* slc) -> float {
    if ( slc->truth.index < 0 ) return -5.0f;

    float count = 0.;
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.pdg == 2212 ) {
        const bool Contained = ( !isnan(prim.end.x) &&
			                           ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
				                          (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
			                           !isnan(prim.end.y) &&
			                           ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                  			         !isnan(prim.end.z) &&
			                           ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );

        if ( Contained ) count+=1.;
      }
    }

    return count;
  });

const Var kProtonMult_True_PostFSI_15MeV([](const caf::SRSliceProxy* slc) -> float {
    if ( slc->truth.index < 0 ) return -5.0f;

    float count = 0.;
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.pdg == 2212 ) {
	      float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
	      //float Et = psquared/(2.*0.938272); // estimate...
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
	      if ( Et > 0.015 ) {
          const bool Contained = ( !isnan(prim.end.x) &&
			                             ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
				                            (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
			                             !isnan(prim.end.y) &&
			                             ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                  			           !isnan(prim.end.z) &&
			                             ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );

          if ( Contained ) count+=1.;
        }
      }
    }

    return count;
  });

const Var kProtonMult_True_PostFSI_100MeV([](const caf::SRSliceProxy* slc) -> float {
    if ( slc->truth.index < 0 ) return -5.0f;

    float count = 0.;
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.pdg == 2212 ) {
        float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
        //float Et = psquared/(2.*0.938272); // estimate...
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
        if ( Et > 0.100 ) {
          const bool Contained = ( !isnan(prim.end.x) &&
			                             ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
				                            (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
			                             !isnan(prim.end.y) &&
			                             ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                  			           !isnan(prim.end.z) &&
			                             ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );

          if ( Contained ) count+=1.;
        }
      }
    }

    return count;
  });

// ---- reco w/o stubs
const Var kProtonMult_Reco([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kPTrackIndNew(slc);
    int countP = 0;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue;  
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      float angle = -5.0;
      if ( primaryInd >= 0 ) {
	const unsigned int idxPrim = (unsigned int)primaryInd;
	TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
	TVector3 pDir( slc->reco.trk[idxTrk].dir.x, slc->reco.trk[idxTrk].dir.y, slc->reco.trk[idxTrk].dir.z );
	angle = TMath::Cos(muDir.Angle(pDir));
      }

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 ) {
	countP += 1;
      }
    }

    return countP;
  });

// ---- reco w/ stubs
const Var kProtonMult_Reco_wStubs([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kPTrackIndNew(slc);
    int countP = 0;

    std::vector< int > proton_pfpids;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; 
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      float angle = -5.0;
      if ( primaryInd >= 0 ) {
        const unsigned int idxPrim = (unsigned int)primaryInd;
        TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
        TVector3 pDir( slc->reco.trk[idxTrk].dir.x, slc->reco.trk[idxTrk].dir.y, slc->reco.trk[idxTrk].dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 ) {
        countP += 1;
	      proton_pfpids.push_back( trk.pfp.id );
      }
    }

    int countPStub = 0;
    for ( unsigned int idxStub = 0; idxStub < slc->reco.stub.size(); ++idxStub ) {
      // Check if stub within 10cm of vertex, contained, and does NOT match one of our already
      // selected proton PFParticles.

      auto const& stub = slc->reco.stub.at(idxStub);
      const float Atslc = std::hypot(slc->vertex.x - stub.vtx.x,
                                     slc->vertex.y - stub.vtx.y,
                                     slc->vertex.z - stub.vtx.z);
      const bool Contained = ( !isnan(stub.end.x) &&
                               ((stub.end.x < -61.94 - 10 && stub.end.x > -358.49 + 10) ||
                                (stub.end.x >  61.94 + 10 && stub.end.x <  358.49 - 10)) &&
                               !isnan(stub.end.y) &&
                               ( stub.end.y > -181.86 + 10 && stub.end.y < 134.96 - 10 ) &&
                               !isnan(stub.end.z) &&
                               ( stub.end.z > -894.95 + 10 && stub.end.z < 894.95 - 10 ) );

      if ( Atslc >= 10.0 && !Contained ) continue;

      // Gray makes the following 2d cut on proton ID for stubs:
      // length <= 0.5cm , dQ/dx > 3.5e5 e/cm
      // else: length <= 2cm, dQ/dx > 3e5 e/cm
      // else: length <= 3cm, dQ/dx > 2.5e5 e/cm
      const float stubLen = sqrt( std::pow(stub.end.x-stub.vtx.x,2) + std::pow(stub.end.y-stub.vtx.y,2) + std::pow(stub.end.z-stub.vtx.z,2) );
      float stubDQ = 0.;
      for ( auto const& plane : stub.planes ) {
        if ( plane.p != 2 ) continue;
        for ( auto const& hit : plane.hits ) stubDQ += hit.charge;
      }
      const float stubDQDx = stubDQ / stubLen;
      if ( stubLen <= 0. ) continue;
      bool passes2Dcut = false;
      if ( stubLen <= 0.5 && stubDQDx > 3.5e5 ) passes2Dcut = true;
      else if ( stubLen <= 2.0 && stubDQDx > 3.0e5 ) passes2Dcut = true;
      else if ( stubLen <= 3.0 && stubDQDx > 2.5e5 ) passes2Dcut = true;
      if ( !passes2Dcut ) continue;
 
      // Now check if it matches an already identified proton PFP based on tracks
      bool matchesProtonPFP = false;
      if ( stub.pfpid >= 0 ) {
        for ( auto const& pfpid : proton_pfpids ) {
          if ( pfpid == stub.pfpid ) matchesProtonPFP = true;
        }
      }
      if ( matchesProtonPFP ) continue;

      // Gray said at the moment the simple selection/analysis of the stubs is not really trustworthy beyond 1 stub... how often would we pick a second?
      countPStub += 1;
      break;
    }

    //std::cout << "countP stub: " << countPStub << std::endl;

    return countP + countPStub;
  });

const SpillVar kTrueProtonStubDQDx( [](const caf::SRSpillProxy *sr) -> double {
    double maxLen = -1.;
    double maxLenDQDx = -1.;

    for ( auto const& slc : sr->slc ) {
      for ( auto const& stub : slc.reco.stub ) {
        if ( stub.truth.p.pdg != 2212 ) continue;
        float p = sqrt(std::pow( stub.truth.p.startp.x, 2 ) + std::pow( stub.truth.p.startp.y, 2 ) + std::pow( stub.truth.p.startp.z, 2 ));
        if ( p > 0.3 ) continue;
        const double stubLen = sqrt( std::pow(stub.end.x-stub.vtx.x,2) + std::pow(stub.end.y-stub.vtx.y,2) + std::pow(stub.end.z-stub.vtx.z,2) );
        double stubDQ = 0.;
        for ( auto const& plane : stub.planes ) {
          if ( plane.p != 2 ) continue;
          for ( auto const& hit : plane.hits ) stubDQ += hit.charge;
        }
        if ( stubDQ < std::numeric_limits<double>::epsilon() ||
             stubLen < std::numeric_limits<double>::epsilon() ) continue;
        const double dqdx = stubDQ / stubLen;
        if ( stubLen > maxLen ) {
          maxLen = stubLen;
          maxLenDQDx = dqdx;
        }
      }
    }

    if ( maxLen < 0. ) return -5.;

    return maxLenDQDx;
  });

const SpillVar kTrueProtonStubLen( [](const caf::SRSpillProxy *sr) -> double {
    double maxLen = -1.;

    for ( auto const& slc : sr->slc ) {
      for ( auto const& stub : slc.reco.stub ) {
        if ( stub.truth.p.pdg != 2212 ) continue;
        float p = sqrt(std::pow( stub.truth.p.startp.x, 2 ) + std::pow( stub.truth.p.startp.y, 2 ) + std::pow( stub.truth.p.startp.z, 2 ));
        if ( p > 0.3 ) continue;
        const double stubLen = sqrt( std::pow(stub.end.x-stub.vtx.x,2) + std::pow(stub.end.y-stub.vtx.y,2) + std::pow(stub.end.z-stub.vtx.z,2) );
        if ( stubLen > maxLen ) maxLen = stubLen;
      }
    }

    if ( maxLen < 0. ) return -5.;

    return maxLen;
  });

const SpillVar kNotProtonStubDQDx( [](const caf::SRSpillProxy *sr) -> double {
    double maxLen = -1.;
    double maxLenDQDx = -1.;

    for ( auto const& slc : sr->slc ) {
      for ( auto const& stub : slc.reco.stub ) {
        if ( stub.truth.p.pdg == 2212 ) continue;
        const double stubLen = sqrt( std::pow(stub.end.x-stub.vtx.x,2) + std::pow(stub.end.y-stub.vtx.y,2) + std::pow(stub.end.z-stub.vtx.z,2) );
        double stubDQ = 0.;
        for ( auto const& plane : stub.planes ) {
          if ( plane.p != 2 ) continue;
          for ( auto const& hit : plane.hits ) stubDQ += hit.charge;
        }
        if ( stubDQ < std::numeric_limits<double>::epsilon() ||
             stubLen < std::numeric_limits<double>::epsilon() ) continue;
        const double dqdx = stubDQ / stubLen;
        if ( stubLen > maxLen ) {
          maxLen = stubLen;
          maxLenDQDx = dqdx;
        }
      }
    }

    if ( maxLen < 0. ) return -5.;

    return maxLenDQDx;
  });

const SpillVar kNotProtonStubLen( [](const caf::SRSpillProxy *sr) -> double {
    double maxLen = -1.;

    for ( auto const& slc : sr->slc ) {
      for ( auto const& stub : slc.reco.stub ) {
        if ( stub.truth.p.pdg == 2212 ) continue;
        const double stubLen = sqrt( std::pow(stub.end.x-stub.vtx.x,2) + std::pow(stub.end.y-stub.vtx.y,2) + std::pow(stub.end.z-stub.vtx.z,2) );
        if ( stubLen > maxLen ) maxLen = stubLen;
      }
    }

    if ( maxLen < 0. ) return -5.;

    return maxLen;
  });


// cheated version with some reco
const Var kProtonMult_Reco_Cheated([](const caf::SRSliceProxy* slc) -> int {
    int countP = 0;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      auto const& trk = slc->reco.trk.at(idxTrk);
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if ( Atslc < 10.0 && slc->reco.trk[idxTrk].pfp.parent_is_primary && slc->reco.trk[idxTrk].truth.p.pdg == 2212  ) countP+=1;
    }

    for ( unsigned int idxShw = 0; idxShw < slc->reco.shw.size(); ++idxShw ) {
      auto const& shw = slc->reco.shw.at(idxShw);
      const float Atslc = std::hypot(slc->vertex.x - shw.start.x,
                                     slc->vertex.y - shw.start.y,
                                     slc->vertex.z - shw.start.z);
      if ( Atslc < 10.0 && slc->reco.shw[idxShw].pfp.parent_is_primary && slc->reco.shw[idxShw].truth.p.pdg == 2212  ) countP+=1;
    }

    return countP;
  });


// Looking at the true protons from neutrinos on a SPILL level
const SpillMultiVar kTrueProtonsKineticEnergy( [](const caf::SRSpillProxy *sr) {
  std::vector<int> g4ids;
  std::vector<double> kes;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 ){
        g4ids.push_back( prim.G4ID );
        float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
        kes.push_back( Et );
      }
    }
  }

  return kes;
});

const SpillMultiVar kTrueProtonsKineticEnergy_CheatedReco( [](const caf::SRSpillProxy *sr) {
  std::vector<double> protonKEs;

  std::vector<int> g4ids;
  std::map<int, double> g4idKEs;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 ){
        g4ids.push_back( prim.G4ID );
        float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
        g4idKEs[prim.G4ID] = Et;
      }
    }
  }

  if ( g4ids.size() == 0 ) return protonKEs;

  std::map< int, unsigned int > g4idFound;
  for ( unsigned int i=0; i<g4ids.size(); ++i ) g4idFound[ g4ids[i] ] = 0;

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
    for ( auto const& shw : slc.reco.shw ) {
      if ( g4idFound.find( shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[shw.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) protonKEs.push_back( g4idKEs[g4id] );
  }

  return protonKEs;
});

const SpillMultiVar kTrueProtonsKineticEnergy_RES( [](const caf::SRSpillProxy *sr) {
  std::vector<int> g4ids;
  std::vector<double> kes;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) ||
         nu.genie_mode!=caf::kRes )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 ){
        g4ids.push_back( prim.G4ID );
        float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
        kes.push_back( Et );
      }
    }
  }

  return kes;
});

const SpillMultiVar kTrueProtonsKineticEnergy_RES_CheatedReco( [](const caf::SRSpillProxy *sr) {
  std::vector<double> protonKEs;

  std::vector<int> g4ids;
  std::map<int, double> g4idKEs;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) ||
         nu.genie_mode!=caf::kRes )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 ){
        g4ids.push_back( prim.G4ID );
        float psquared = std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 );
        float Et = sqrt( psquared + std::pow(0.938272,2) ) - 0.938272;
        g4idKEs[prim.G4ID] = Et;
      }
    }
  }

  if ( g4ids.size() == 0 ) return protonKEs;

  std::map< int, unsigned int > g4idFound;
  for ( unsigned int i=0; i<g4ids.size(); ++i ) g4idFound[ g4ids[i] ] = 0;

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
    for ( auto const& shw : slc.reco.shw ) {
      if ( g4idFound.find( shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[shw.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) protonKEs.push_back( g4idKEs[g4id] );
  }

  return protonKEs;
});

// BH: didn't add true containment to this check.. Let's add that for the momentum version...
// Same but with momentum not kinetic energy
const SpillMultiVar kTrueProtonsMomentum( [](const caf::SRSpillProxy *sr) {
  std::vector<int> g4ids;
  std::vector<double> ps;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        g4ids.push_back( prim.G4ID );
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        ps.push_back( p );
      }
    }
  }

  return ps;
});

const SpillMultiVar kTrueProtonsMomentum_CheatedReco( [](const caf::SRSpillProxy *sr) {
  std::vector<double> ps;

  std::vector<int> g4ids;
  std::map<int, double> g4idPs;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        g4ids.push_back( prim.G4ID );
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        g4idPs[prim.G4ID] = p;
      }
    }
  }

  if ( g4ids.size() == 0 ) return ps;

  std::map< int, unsigned int > g4idFound;
  for ( unsigned int i=0; i<g4ids.size(); ++i ) g4idFound[ g4ids[i] ] = 0;

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
    for ( auto const& shw : slc.reco.shw ) {
      if ( g4idFound.find( shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[shw.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idPs[g4id] );
  }

  return ps;
});

const SpillMultiVar kTrueProtonsMomentum_CheatedRecoTrk( [](const caf::SRSpillProxy *sr) {
  std::vector<double> ps;

  std::vector<int> g4ids;
  std::map<int, double> g4idPs;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        g4ids.push_back( prim.G4ID );
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        g4idPs[prim.G4ID] = p;
      }
    }
  }

  if ( g4ids.size() == 0 ) return ps;

  std::map< int, unsigned int > g4idFound;
  for ( unsigned int i=0; i<g4ids.size(); ++i ) g4idFound[ g4ids[i] ] = 0;

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idPs[g4id] );
  }

  return ps;
});

// Not cheated reco -- applies track containment and chi2 cut...
const SpillMultiVar kTrueProtonsMomentum_Tracked( [](const caf::SRSpillProxy *sr) {
  std::vector<double> ps;

  std::vector<int> g4ids;
  std::map<int, double> g4idPs;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        g4ids.push_back( prim.G4ID );
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        g4idPs[prim.G4ID] = p;
      }
    }
  }

  if ( g4ids.size() == 0 ) return ps;

  std::map< int, unsigned int > g4idFound;
  for ( unsigned int i=0; i<g4ids.size(); ++i ) g4idFound[ g4ids[i] ] = 0;

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( trk.bestplane == -1 ) continue;
      // Check cuts
      const bool Contained = ( !isnan(trk.end.x) &&
			                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
			                        	(trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			                         !isnan(trk.end.y) &&
			                         ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			                         !isnan(trk.end.z) &&
			                         ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( !Contained || Chi2Proton > 100. || Chi2Muon < 30. ) continue;
      // const float pida = trk.chi2pid[trk.bestplane].pida;
      //if ( Contained && ( pida > 15 || ( Chi2Proton <= 100. && Chi2Muon >= 30. ) ) ) {
      // See if it matches our proton G4ID
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
      //}
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idPs[g4id] );
  }

  return ps;
});

const SpillMultiVar kTrueProtonsRTLen_CheatedRecoTrk( [](const caf::SRSpillProxy *sr) {
  std::vector<double> ps;

  std::vector<int> g4ids;
  std::map<int, double> g4idLens;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        g4ids.push_back( prim.G4ID );
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        g4idLens[prim.G4ID] = prim.length;
      }
    }
  }

  if ( g4ids.size() == 0 ) return ps;

  std::map< int, unsigned int > g4idFound;
  std::map< int, float > g4idSelTrkLen;
  for ( unsigned int i=0; i<g4ids.size(); ++i ){
    g4idFound[ g4ids[i] ] = 0;
    g4idSelTrkLen[ g4ids[i] ] = 0.;
  }

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ){
        g4idFound[trk.truth.p.G4ID]+=1;
        if ( trk.len > g4idSelTrkLen[ trk.truth.p.G4ID ] ) g4idSelTrkLen[ trk.truth.p.G4ID ] = trk.len;
      }
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idSelTrkLen[g4id] / g4idLens[g4id] );
  }

  return ps;
});

const SpillMultiVar kTrueProtonsRTLen_CheatedRecoTrk_HighE( [](const caf::SRSpillProxy *sr) {
  std::vector<double> ps;

  std::vector<int> g4ids;
  std::map<int, double> g4idLens;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    for ( auto const& prim : nu.prim ) {
      if ( prim.pdg == 2212 && prim.contained ){
        // What if we define containment with 10cm like in reco?
        const bool TContained = (!isnan(prim.end.x) &&
                                 ((prim.end.x < -61.94 - 10 && prim.end.x > -358.49 + 10) ||
                                  (prim.end.x >  61.94 + 10 && prim.end.x <  358.49 - 10)) &&
                                 !isnan(prim.end.y) &&
                                 ( prim.end.y > -181.86 + 10 && prim.end.y < 134.96 - 10 ) &&
                                 !isnan(prim.end.z) &&
                                 ( prim.end.z > -894.95 + 10 && prim.end.z < 894.95 - 10 ) );
        if ( !TContained ) continue;
        // -----------------------------------------------------
        float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
        if ( p < 1.0 ) continue; // minimum |p| = 1 GeV/c here...
        g4ids.push_back( prim.G4ID );
        g4idLens[prim.G4ID] = prim.length;
      }
    }
  }

  if ( g4ids.size() == 0 ) return ps;

  std::map< int, unsigned int > g4idFound;
  std::map< int, float > g4idSelTrkLen;
  for ( unsigned int i=0; i<g4ids.size(); ++i ){
    g4idFound[ g4ids[i] ] = 0;
    g4idSelTrkLen[ g4ids[i] ] = 0.;
  }

  for ( auto const& slc : sr->slc ) {
    for ( auto const& trk : slc.reco.trk ) {
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ){
        g4idFound[trk.truth.p.G4ID]+=1;
        if ( trk.len > g4idSelTrkLen[ trk.truth.p.G4ID ] ) g4idSelTrkLen[ trk.truth.p.G4ID ] = trk.len;
      }
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idSelTrkLen[g4id] / g4idLens[g4id] );
  }

  return ps;
});





// The base cut in the committed helper would be
//
//   const Cut kNuMuCC_FullSelection = ( kCryo0 && kRFiducial && kNotClearCosmic && kNuScore && kFMScore && kPTrack );
//
// and we want to
//    * ditch the kCryo0
//    * use our updated FV cut on the vertex here
//    * keep kNotClearCosmic
//    * REPLACE NuScore cut with a cut on the Longest Track Y direction in CR hypothesis --> make that cut here...
//    * use FMScore but change the threshold to 12
//    * ADD an FMTime cut...
//    * Use kPTrack BUT a version that checks against the new values...
//    * ADD check for a secondary track that is within 10cm of vertex and contained -- do NOT check its chi2_proton

// Copying the one cut that isn't altered, so I can get rid of the includes...
const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });
// ---------------

const Cut kPTrackNew([](const caf::SRSliceProxy* slc) {
    return ( kPTrackIndNew(slc) >= 0 );
  });

const Cut kRFiducialNew([](const caf::SRSliceProxy* sr)
			{
			  const auto& vtx = sr->vertex;

			  if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

			  return (( ( vtx.x < -61.94 - 25 && vtx.x > -358.49 + 25 ) ||
				    ( vtx.x >  61.94 + 25 && vtx.x <  358.49 - 25 )) &&
				  (( vtx.y > -181.86 + 25 && vtx.y < 134.96 - 25 ) &&
				   ( vtx.z > -894.95 + 30 && vtx.z < 894.95 - 50 ) ));
			});

// Typically here copied from Numu Vars and/or Cuts and altered
const Cut kFMScoreNuMI([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 );
  });

const Cut kFMTimeNuMI([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.time) && slc->fmatch.time > -0.2 && slc->fmatch.time < 9.9 );
  });

const Cut kCutCRLongTrkDirY([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->nuid.crlongtrkdiry) && slc->nuid.crlongtrkdiry > -0.7 );
  });

const Cut kProtonTrack([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
				                             slc->vertex.y - trk.start.y,
				                             slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
			                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
			                        	(trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			                         !isnan(trk.end.y) &&
			                         ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			                         !isnan(trk.end.z) &&
			                         ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( Chi2Proton > 100. || Chi2Muon < 30. ) continue;
      ////////////

      TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
      TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && trk.len > maxLength && TMath::Cos(muDir.Angle(pDir)) >= -0.9 ) {
      	maxLength = trk.len;
      	idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return false;

    return true;
  });

const Cut kProtonTrack400MeV([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
				                             slc->vertex.y - trk.start.y,
				                             slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
			                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
			                        	(trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			                         !isnan(trk.end.y) &&
			                         ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			                         !isnan(trk.end.z) &&
			                         ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( Chi2Proton > 100. || Chi2Muon < 30. ) continue;
      ////////////

      TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
      TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );

      float protonp = -1.;
      if(Contained) protonp = trk.rangeP.p_proton;
      //else {
  	  //  protonp = trk.mcsP.fwdP_proton;
	    //}

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && trk.len > maxLength && TMath::Cos(muDir.Angle(pDir)) >= -0.9 && protonp > 0.4 ) {
      	maxLength = trk.len;
      	idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return false;

    return true;
  });

//const Cut kNuMISelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kFMScoreNuMI && kFMTimeNuMI && kPTrackNew && kProtonTrack );
// TPC-only
const Cut kNuMISelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kProtonTrack );
const Cut kNuMISelectionContained = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kProtonTrack && kRecoMuonContained );
const Cut kNuMISelection_TrueSignal = kCutTrueSig && kNuMISelection;

const Cut kNuMISelection_NoTrkDirY = ( kRFiducialNew && kNotClearCosmic && kPTrackNew && kProtonTrack );

// 1Mu1P (well really 1MuNP with a minimum P threshold) selection
const Cut kNuMI1Mu1PSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kProtonTrack400MeV );

// Version without proton cut for investigation on the proton cut...
const Cut kNuMISelectionWOProton = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew );

/////// Experiment with a COHERENT pion production sample
// Make a pion cut similar to the proton cut above, but swap out the chi2proton -> chi2pion and chi2muon -> chi2proton. And require == 2 trks within 10cm of vtx
const Cut kOnlyScdyPionTrack([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    unsigned int primTrks = 1; // start at 1 because muon

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Pion = trk.chi2pid[trk.bestplane].chi2_pion;
      if ( Chi2Pion > 100. || Chi2Proton < 30. ) continue;
      ////////////

      if ( Atslc < 10.0 ) primTrks+=1;

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( primTrks != 2 ) return false;
    if ( maxLength < 0. ) return false;

    TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
    TVector3 piDir( slc->reco.trk[idxScdy].dir.x, slc->reco.trk[idxScdy].dir.y, slc->reco.trk[idxScdy].dir.z );
    if ( TMath::Cos(muDir.Angle(piDir)) < -0.9 ) return false;
    return true;
  }); // kOnlyScdyPionTrack

// No Showers Cut:
const Cut kNoAppreciableShower([](const caf::SRSliceProxy* slc) {
    unsigned int appreciableShowers = 0;
    for ( auto const& shw : slc->reco.shw ) 
    {
      const float Atslc = std::hypot(slc->vertex.x - shw.start.x,
                                     slc->vertex.y - shw.start.y,
                                     slc->vertex.z - shw.start.z);

      if ( Atslc < 15. && shw.bestplane_energy > 0.4 ) appreciableShowers+=1;
    }

    if ( appreciableShowers == 0 ) return true;
    return false;
  });

// No protons
const Cut kNoProtonTrack([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( Chi2Proton > 100. || Chi2Muon < 30. ) continue;
      ////////////

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return true;
    return false;
  });

const Cut kCheatedNoProton([](const caf::SRSliceProxy* slc) {
    if (slc->truth.index < 0) return false; // not a Nu
    for ( auto const& part : slc->truth.prim ) {
      if ( part.pdg == 2212 && sqrt( std::pow(part.startp.x,2) + 
				     std::pow(part.startp.y,2) + 
				     std::pow(part.startp.z,2) ) > .050 ) return false;
    }

    return true;
  });

const Cut kCheatedNoProtonLowE([](const caf::SRSliceProxy* slc) {
    if (slc->truth.index < 0) return false; // not a Nu
    for ( auto const& part : slc->truth.prim ) {
      if ( part.pdg == 2212 && sqrt( std::pow(part.startp.x,2) +
                                     std::pow(part.startp.y,2) +
                                     std::pow(part.startp.z,2) ) > .015 ) return false;
    }

    return true;
  });

const Var kTCoherent([](const caf::SRSliceProxy* slc) -> double {
    // See https://arxiv.org/pdf/1409.3835.pdf

    // Assuming we've already made the cut to things with 2 primary tracks and have 2 contained tracks,
    // then we can do things a bit more simply here...
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return -5.0;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue;

      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      if ( Atslc < 10.0 ) {
	idxScdy = idxTrk;
	break;
      }
    }

    // Muon:
    double pmuon = slc->reco.trk[idxPrim].rangeP.p_muon;
    double Emuon = sqrt((.105658*.105658) + (pmuon*pmuon));
    TVector3 vecpmuon(pmuon*slc->reco.trk[idxPrim].dir.x,pmuon*slc->reco.trk[idxPrim].dir.y,pmuon*slc->reco.trk[idxPrim].dir.z);

    double plmuon = vecpmuon.Dot(NuDirection_NuMI); // longitundinal momentum
    TVector3 ptmuon = vecpmuon.Cross(NuDirection_NuMI); // transverse momentum

    // Pion:
    double ppion = slc->reco.trk[idxScdy].rangeP.p_pion;
    double Epion = sqrt((.139571*.139571) + (ppion*ppion));
    TVector3 vecppion(ppion*slc->reco.trk[idxScdy].dir.x,ppion*slc->reco.trk[idxScdy].dir.y,ppion*slc->reco.trk[idxScdy].dir.z);

    double plpion = vecppion.Dot(NuDirection_NuMI);
    TVector3 ptpion = vecppion.Cross(NuDirection_NuMI);

    // |t| calculation:
    double term1 = std::pow((Emuon-plmuon)+(Epion-plpion),2);
    ptmuon+=ptpion;
    double term2 = ptmuon.Dot(ptmuon);
    return term1+term2;
  });

// Swap the proton cut for a pion cut, and force the ! of the proton cut
//const Cut kNuMICOHSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kOnlyScdyPionTrack && !kProtonTrack );
//const Cut kNuMICOHSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kOnlyScdyPionTrack && kNoAppreciableShower && kRecoMuonContained && !kProtonTrack );
const Cut kNuMICOHSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kOnlyScdyPionTrack && kNoAppreciableShower && kRecoMuonContained && kNoProtonTrack );

const Cut kNuMICOHSelection_CheatedNoP = kNuMICOHSelection && kCheatedNoProton;
const Cut kNuMICOHSelection_CheatedNoPLowE = kNuMICOHSelection && kCheatedNoProtonLowE;



// Modify cut looking for ONLY a secondary pion to instead be at LEAST a secondary pion
const Cut kHasScdyPionTrack([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Pion = trk.chi2pid[trk.bestplane].chi2_pion;
      if ( Chi2Pion > 100. || Chi2Proton < 30. ) continue;
      ////////////

      if ( Atslc < 10.0 && trk.pfp.parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return false;

    TVector3 muDir( slc->reco.trk[idxPrim].dir.x, slc->reco.trk[idxPrim].dir.y, slc->reco.trk[idxPrim].dir.z );
    TVector3 piDir( slc->reco.trk[idxScdy].dir.x, slc->reco.trk[idxScdy].dir.y, slc->reco.trk[idxScdy].dir.z );
    if ( TMath::Cos(muDir.Angle(piDir)) < -0.9 ) return false;
    return true;
  }); // kHasScdyPionTrack

const Cut kThreePrimaryTracks([](const caf::SRSliceProxy* slc) {
    unsigned int nPrim=0;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      auto const& trk = slc->reco.trk.at(idxTrk);
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if ( Atslc < 10.0 && trk.pfp.parent_is_primary ) nPrim+=1;
    }

    return nPrim==3;
  });

const Cut kNuMIRESSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kHasScdyPionTrack && kNoAppreciableShower && kProtonTrack && kThreePrimaryTracks );

/////// Gray's proposed sample
const Cut kGraysProposedSampleCut([](const caf::SRSliceProxy* slc) {
    if (slc->is_clear_cosmic) return false;
    if (slc->reco.trk.size() < 3) return false;

    unsigned int nPrimTracks = 1; // think after talking to Jaesung this needs to be 1 - we found this bug while comparing cuts
    bool oneIsProton = false;

    for ( auto const& trk : slc->reco.trk ) {
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      if (Atslc < 10.) nPrimTracks+=1;
      else continue;

      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      if ( Contained && Chi2Proton < 60. && Chi2Muon > 40. ) oneIsProton = true;
    }

    if ( oneIsProton && nPrimTracks >= 3 ) return true;
    return false;
  });

const Cut kGraysProposedSampleCutWithMuonTrack([](const caf::SRSliceProxy* slc) {
    if (slc->is_clear_cosmic) return false;
    if (slc->reco.trk.size() < 3) return false;

    unsigned int nPrimTracks = 1; // think after talking to Jaesung this needs to be 1 - we found this bug while comparing cuts
    bool oneIsProton = false;

    int idxMuon = kPTrackIndNew(slc);
    if ( idxMuon < 0 ) return false;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == idxMuon ) continue; // need a different track...
  
      auto const& trk = slc->reco.trk.at(idxTrk);
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      if (Atslc < 10.) nPrimTracks+=1;
      else continue;

      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      if ( Contained && Chi2Proton < 60. && Chi2Muon > 40. ) oneIsProton = true;
    }

    if ( oneIsProton && nPrimTracks >= 3 ) return true;
    return false;
  });

const Var kProtonTrackIndGraysSample([](const caf::SRSliceProxy* slc) -> int {
    int index = -1;
    double longestCandidate = 0.;

    int idxMuon = kPTrackIndNew(slc);
    if ( idxMuon < 0 ) return -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( (int)idxTrk == idxMuon ) continue; // need a different track...

      auto const& trk = slc->reco.trk.at(idxTrk);
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      if (Atslc >= 10.) continue;

      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      if ( Contained && Chi2Proton < 60. && Chi2Muon > 40. && trk.len > longestCandidate ) {
        index = idxTrk;
        longestCandidate = trk.len;
      }
    }

    return index;
  });

const Var kRecoProtonPGraysSample([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);

  if ( kProtonTrackIndGraysSample(slc) >= 0 )
  {
    auto const& trk = slc->reco.trk.at(kProtonTrackIndGraysSample(slc));
    const bool Contained = ( !isnan(trk.end.x) &&
                             ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                              (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                             !isnan(trk.end.y) &&
                             ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                             !isnan(trk.end.z) &&
                             ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
    if(Contained) p = trk.rangeP.p_proton;
    else {
	    std::cout << "Currently kProtonTrackIndGraysSample requires a contained proton... Why am I trying to use MCS here??" << std::endl;
	    p = trk.mcsP.fwdP_proton;
	  }
  }

  return p;
});

const Cut kRecoProtonIsTrueProtonGraysSample([](const caf::SRSliceProxy* slc) {
  if ( kProtonTrackIndGraysSample(slc) >= 0 )
  {
    auto const& trk = slc->reco.trk.at(kProtonTrackIndGraysSample(slc));

    if ( trk.truth.p.pdg == 2212 ) return true;
    else return false;
  }

  return false;
});


//////////////////////////////////////////////////////////////////////
// Spill MultiVars for true signal E_nu and for selected slices...
//////////////////////////////////////////////////////////////////////

const SpillMultiVar kTrueSignalEnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalEnu;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          break;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          break;
        }
      }
    }

		if ( signalLepton ) signalEnu.push_back( nu.E );
  }

  return signalEnu;
});

const SpillMultiVar kTrueSignalwProtonMuMom ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalMuMom;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;
    bool signalHadron = false;

    float muonlen = -1.;
    float muonp = -1.;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. && prim.length > muonlen ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained && prim.length > muonlen ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 2212 && !signalHadron ) {
          //if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 && p < 1.0 )
            signalHadron = true;
          continue;
        }
      }
    }

		if ( signalLepton && signalHadron && muonlen > 0. ) signalMuMom.push_back( muonp );
  }

  return signalMuMom;
});

const SpillMultiVar kTrueSignalLnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalLnu;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          break;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          break;
        }
      }
    }

		if ( signalLepton ) signalLnu.push_back( nu.baseline );
  }

  return signalLnu;
});

const SpillMultiVar kTrueSignalLOverEnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalLOverEnu;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          break;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          break;
        }
      }
    }

		if ( signalLepton ) signalLOverEnu.push_back( nu.baseline / (1000.*nu.E) );
  }

  return signalLOverEnu;
});

// Applying selections

const SpillMultiVar kTrueSignalEnuSelected ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalEnu;
  std::map<int, double> signalIndexE;
  unsigned int count = 0;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          break;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          break;
        }
      }
    }

		if ( signalLepton ) {
      signalIndexE[ nu.index ] = nu.E;
      count+=1;
    }
  }

  if ( count == 0 ) return signalEnu;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( signalIndexE.find( slc.truth.index ) == signalIndexE.end() ) continue;

    // Check if slice passes cuts:
    if ( kRFiducialNew(&slc) &&
         kNotClearCosmic(&slc) &&
         kCutCRLongTrkDirY(&slc) &&
         kPTrackNew(&slc) &&
         kProtonTrack(&slc) ) {
      signalEnu.push_back( slc.truth.E );
      signalIndexE.erase( slc.truth.index ); // should prevent duplicates
    }
  }

  return signalEnu;
});

const SpillMultiVar kTrueSignalwProtonMuMomSelected ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalMuMom;
  std::map<int, double> signalIndexMuMom;
  unsigned int count = 0;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;
    bool signalHadron = false;

    float muonlen = -1.;
    float muonp = -1.;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 2212 && !signalHadron ) {
          //if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 && p < 1.0 )
            signalHadron = true;
          continue;
        }
      }
    }

		if ( signalLepton && signalHadron && muonlen > 0. ) {
      signalIndexMuMom[ nu.index ] = muonp;
      count+=1;
    }
  }

  if ( count == 0 ) return signalMuMom;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( signalIndexMuMom.find( slc.truth.index ) == signalIndexMuMom.end() ) continue;

    // Check if slice passes cuts:
    if ( kRFiducialNew(&slc) &&
         kNotClearCosmic(&slc) &&
         kCutCRLongTrkDirY(&slc) &&
         kPTrackNew(&slc) &&
         kProtonTrack400MeV(&slc) &&
         kRecoMuonContained(&slc) ) {
      signalMuMom.push_back( signalIndexMuMom.at(slc.truth.index) );
      signalIndexMuMom.erase( slc.truth.index ); // should prevent duplicates
    }
  }

  return signalMuMom;
});

// --- RECO VERSIONS FOR PURITY
const SpillMultiVar kTrueSignalwProtonRecoMuMomSelected ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalMuMom;
  std::map<int, double> signalIndexMuMom;
  unsigned int count = 0;

  for ( auto const& nu : sr->mc.nu ) {
    bool signalLepton = false;
    bool signalHadron = false;

    float muonlen = -1.;
    float muonp = -1.;

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      for ( auto const& prim : nu.prim ) {
        if ( prim.pdg == 13 && prim.length > 100. ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 13 && prim.length > 50. && prim.contained ) {
          signalLepton = true;
          muonlen = prim.length;
          muonp = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          continue;
        }
        else if ( prim.pdg == 2212 && !signalHadron ) {
          //if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 && p < 1.0 )
            signalHadron = true;
          continue;
        }
      }
    }

		if ( signalLepton && signalHadron && muonlen > 0. ) {
      signalIndexMuMom[ nu.index ] = muonp;
      count+=1;
    }
  }

  if ( count == 0 ) return signalMuMom;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( signalIndexMuMom.find( slc.truth.index ) == signalIndexMuMom.end() ) continue;

    // Check if slice passes cuts:
    if ( kRFiducialNew(&slc) &&
         kNotClearCosmic(&slc) &&
         kCutCRLongTrkDirY(&slc) &&
         kPTrackNew(&slc) &&
         kProtonTrack400MeV(&slc) &&
         kRecoMuonContained(&slc) ) {
      signalMuMom.push_back( kRecoMuonPNew(&slc) );
      signalIndexMuMom.erase( slc.truth.index ); // should prevent duplicates
    }
  }

  return signalMuMom;
});

const SpillMultiVar kAllwProtonRecoMuMomSelected ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalMuMom;

  for ( auto const& slc : sr->slc ) {
    // Check if slice passes cuts:
    if ( kRFiducialNew(&slc) &&
         kNotClearCosmic(&slc) &&
         kCutCRLongTrkDirY(&slc) &&
         kPTrackNew(&slc) &&
         kProtonTrack400MeV(&slc) &&
         kRecoMuonContained(&slc) ) {
      signalMuMom.push_back( kRecoMuonPNew(&slc) );
    }
  }

  return signalMuMom;
});


// Version with all fiducial numu cc
const SpillMultiVar kTrueNumuCCEnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalEnu;

  for ( auto const& nu : sr->mc.nu ) {

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      signalEnu.push_back( nu.E );
    }

  }

  return signalEnu;
});

const SpillMultiVar kTrueNumuCCLnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalLnu;

  for ( auto const& nu : sr->mc.nu ) {

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      signalLnu.push_back( nu.baseline );
    }

  }

  return signalLnu;
});

const SpillMultiVar kTrueNumuCCLOverEnu ( [](const caf::SRSpillProxy *sr) {
  std::vector<double> signalLOverEnu;

  for ( auto const& nu : sr->mc.nu ) {

    if ( abs(nu.pdg) == 14 &&
				 nu.iscc &&
				 !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
				 isInFV(nu.position.x,nu.position.y,nu.position.z) ) {
      signalLOverEnu.push_back( nu.baseline / (1000.*nu.E) );
    }

  }

  return signalLOverEnu;
});







// SpillVars to calcualte ThXZ, ThNuMI, etc.
const SpillVar kSingleMatchCosThXZ ([](const caf::SRSpillProxy* sr) -> double
                                 {
                                   double thXZ = 0.;
                                   double length = 0.;

                                   int g4id = -1;
                                   int cryo = -1;

                                   std::pair<unsigned int, unsigned int> slcTrk;
                                   double maxEnergy = -1.;

                                   for ( auto const& iNu : sr->mc.nu ){
                                     for ( auto const& ipart : iNu.prim ){
                                       if ( ipart.pdg != 13 ) continue;
                                       if ( ipart.start_process!=0 ||
                                            std::fabs(ipart.start.x)<70. || std::fabs(ipart.start.x)>350.||
                                            ipart.start.z < -870. || ipart.start.z > 870. || ipart.start.y < -170. || ipart.start.y > 120. )
                                         continue;
                                       if ( ipart.startT < 0. || ipart.startT >= 9.5 ) continue;
                                       if ( ipart.length > length && ipart.length > 50. ) {
                                         g4id = ipart.G4ID;
                                         length = ipart.length;
                                         cryo = ipart.start.x < 0 ? 0 : 1;
                                          
                                         TVector3 dir(ipart.end.x-ipart.start.x,
                                                      ipart.end.y-ipart.start.y,
                                                      ipart.end.z-ipart.start.z);
                                         dir = dir.Unit();
                                         thXZ = TMath::ATan2(dir.X(),dir.Z());
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thXZ);
                                   //return thXZ;
                                 });

const SpillVar kSingleMatchCosThNuMI ([](const caf::SRSpillProxy* sr) -> double
                                 {
                                   double thNuMI = 0.;
                                   double length = 0.;

                                   int g4id = -1;
                                   int cryo = -1;

                                   std::pair<unsigned int, unsigned int> slcTrk;
                                   double maxEnergy = -1.;

                                   for ( auto const& iNu : sr->mc.nu ){
                                     for ( auto const& ipart : iNu.prim ){
                                       if ( ipart.pdg != 13 ) continue;
                                       if ( ipart.start_process!=0 ||
                                            std::fabs(ipart.start.x)<70. || std::fabs(ipart.start.x)>350.||
                                            ipart.start.z < -870. || ipart.start.z > 870. || ipart.start.y < -170. || ipart.start.y > 120. )
                                         continue;
                                       if ( ipart.startT < 0. || ipart.startT >= 9.5 ) continue;
                                       if ( ipart.length > length && ipart.length > 50. ) {
                                         g4id = ipart.G4ID;
                                         length = ipart.length;
                                         cryo = ipart.start.x < 0 ? 0 : 1;
                                          
                                         TVector3 dir(ipart.end.x-ipart.start.x,
                                                      ipart.end.y-ipart.start.y,
                                                      ipart.end.z-ipart.start.z);
                                         dir = dir.Unit();
                                         thNuMI = dir.Angle(NuDirection_NuMI);
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thNuMI);
                                 });

const SpillVar kSingleMatchCosThNuMINotIsoch ([](const caf::SRSpillProxy* sr) -> double
                                 {
                                   double thNuMI = 0.;
                                   double length = 0.;

                                   int g4id = -1;
                                   int cryo = -1;

                                   std::pair<unsigned int, unsigned int> slcTrk;
                                   double maxEnergy = -1.;

                                   for ( auto const& iNu : sr->mc.nu ){
                                     for ( auto const& ipart : iNu.prim ){
                                       if ( ipart.pdg != 13 ) continue;
                                       if ( ipart.start_process!=0 ||
                                            std::fabs(ipart.start.x)<70. || std::fabs(ipart.start.x)>350.||
                                            ipart.start.z < -870. || ipart.start.z > 870. || ipart.start.y < -170. || ipart.start.y > 120. )
                                         continue;
                                       if ( ipart.startT < 0. || ipart.startT >= 9.5 ) continue;
                                       if ( ipart.length > length && ipart.length > 50. ) {
                                         g4id = ipart.G4ID;
                                         length = ipart.length;
                                         cryo = ipart.start.x < 0 ? 0 : 1;
                                          
                                         TVector3 dir(ipart.end.x-ipart.start.x,
                                                      ipart.end.y-ipart.start.y,
                                                      ipart.end.z-ipart.start.z);
                                         dir = dir.Unit();
                                         double thXZ = (180./TMath::Pi())*TMath::ATan2(dir.X(),dir.Z());
                                         if ( fabs(thXZ) <= 10. ) continue; // very approximate
                                         thNuMI = dir.Angle(NuDirection_NuMI);
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thNuMI);
                                 });

// TODO: NEEDS FIX -- I think here I want like >= 75 and <= 105 or something like this?
const SpillVar kSingleMatchCosThNuMINotIsochNotPerp ([](const caf::SRSpillProxy* sr) -> double
                                 {
                                   double thNuMI = 0.;
                                   double length = 0.;

                                   int g4id = -1;
                                   int cryo = -1;

                                   std::pair<unsigned int, unsigned int> slcTrk;
                                   double maxEnergy = -1.;

                                   for ( auto const& iNu : sr->mc.nu ){
                                     for ( auto const& ipart : iNu.prim ){
                                       if ( ipart.pdg != 13 ) continue;
                                       if ( ipart.start_process!=0 ||
                                            std::fabs(ipart.start.x)<70. || std::fabs(ipart.start.x)>350.||
                                            ipart.start.z < -870. || ipart.start.z > 870. || ipart.start.y < -170. || ipart.start.y > 120. )
                                         continue;
                                       if ( ipart.startT < 0. || ipart.startT >= 9.5 ) continue;
                                       if ( ipart.length > length && ipart.length > 50. ) {
                                         g4id = ipart.G4ID;
                                         length = ipart.length;
                                         cryo = ipart.start.x < 0 ? 0 : 1;
                                          
                                         TVector3 dir(ipart.end.x-ipart.start.x,
                                                      ipart.end.y-ipart.start.y,
                                                      ipart.end.z-ipart.start.z);
                                         dir = dir.Unit();
                                         double thXZ = (180./TMath::Pi())*TMath::ATan2(dir.X(),dir.Z());
                                         if ( fabs(thXZ) <= 10. ) continue; // very approximate
                                         if ( fabs(thXZ) >= 75. ) continue; // very approximate
                                         thNuMI = dir.Angle(NuDirection_NuMI);
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thNuMI);
                                 });


// Attempt to make a 0pion cut

// First attempt bakes it into the Muon selection
const Var kPTrackIndNew_0Pi([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);

    unsigned int nCands(0);

    for (std::size_t i(0); i < slc->reco.trk.size(); ++i)
    {
      auto const& trk = slc->reco.trk.at(i);
      if(trk.bestplane == -1) continue;

      // First we calculate the distance of each track to the slice vertex.
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                      slc->vertex.y - trk.start.y,
                                      slc->vertex.z - trk.start.z);

      // We require that the distance of the track from the slice is less than
      // 10 cm and that the parent of the track has been marked as the primary.
      const bool AtSlice = ( Atslc < 10.0 && trk.pfp.parent_is_primary);

      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      const bool Contained = (!isnan(trk.end.x) &&
                              ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                              !isnan(trk.end.y) &&
                              ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                              !isnan(trk.end.z) &&
                              ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
      const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
      if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) )
        nCands+=1;
      if ( trk.len > Longest )
      {
        Longest = trk.len;
        PTrackInd = i;
      }
    }

    if ( nCands > 1 ) return -1;
    return PTrackInd;
  });

  const Cut kPTrackNew_0Pi([](const caf::SRSliceProxy* slc) {
    return ( kPTrackIndNew_0Pi(slc) >= 0 );
  });

// Second attempt:
// -- No Appreciable shower
// -- Any track that isn't our primary muon or proton needs to be consistent with a proton
//     -- What if we require ONLY two tracks?

const Cut kNoAppreciableShower_0Pi([](const caf::SRSliceProxy* slc) {
    unsigned int appreciableShowers = 0;
    for ( auto const& shw : slc->reco.shw ) 
    {
      const float Atslc = std::hypot(slc->vertex.x - shw.start.x,
                                     slc->vertex.y - shw.start.y,
                                     slc->vertex.z - shw.start.z);

      if ( Atslc < 15. && shw.bestplane_energy > 0.030 ) appreciableShowers+=1;
    }

    if ( appreciableShowers == 0 ) return true;
    return false;
  });

const Cut kCutScdyAndNonProtonShowers_0Pi([](const caf::SRSliceProxy* slc) {
    unsigned int nonPrimaryShw = 0;
    unsigned int nonProtonLikeShowers = 0;

    for ( auto const& shw : slc->reco.shw ) 
    {
      const float Atslc = std::hypot(slc->vertex.x - shw.start.x,
                                     slc->vertex.y - shw.start.y,
                                     slc->vertex.z - shw.start.z);

      if ( Atslc < 15. && shw.bestplane_energy > 0.030 ) {
        if ( !shw.pfp.parent_is_primary ) nonPrimaryShw+=1;
        else if ( shw.density < 7.0 ) nonProtonLikeShowers+=1;
      }
    }

    if ( nonPrimaryShw==0 && nonProtonLikeShowers==0 ) return true;
    return false;
  });

const Cut kNonMuonTracksAreProtons([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxPrim = (unsigned int)primaryInd;

    // Check to see if any of the non-muon candidates is inconsistent with a proton 
    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // muon should be muon-like...

      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue; // can't determine chi2 for this track... sigh.

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
				                             slc->vertex.y - trk.start.y,
				                             slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
			                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
			                        	(trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			                         !isnan(trk.end.y) &&
			                         ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			                         !isnan(trk.end.z) &&
			                         ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      if ( Atslc >= 10.0 || !trk.pfp.parent_is_primary ) continue; // not a primary so we choose to move on

      // Now check if it's actually proton like and return false if not...
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( Chi2Proton > 100. || Chi2Muon < 30. ) return false;
    }

    // All our primary tracks are proton-like, we're good!
    return true;
  });

// const Cut kNuMISelection_0Pi = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew_0Pi && kProtonTrack );

// This looked good but I am doubtful of the no appreciable showers cut...
// const Cut kNuMISelection_0Pi = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && 
//                                  kPTrackNew && kProtonTrack && kNoAppreciableShower_0Pi && kNonMuonTracksAreProtons );

// const Cut kNuMISelection_0Pi = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && 
//                                  kPTrackNew && kProtonTrack && kNonMuonTracksAreProtons );

const Cut kNuMISelection_0Pi = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && 
                                 kPTrackNew && kProtonTrack && kCutScdyAndNonProtonShowers_0Pi && kNonMuonTracksAreProtons );

const Cut kNuMISelection_0Pi_NoTrkDirY = ( kRFiducialNew && kNotClearCosmic && kPTrackNew && kProtonTrack &&
                                           kCutScdyAndNonProtonShowers_0Pi && kNonMuonTracksAreProtons );

/////////
// EXPLORATION

const Cut kNonMuonTracksAreProtons__PRINT([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    int secondaryInd = kPTrackIndNewProton(slc);
    if ( secondaryInd < 0 ) return false;
    unsigned int idxScdy = (unsigned int)secondaryInd;

    bool shouldCutEvt = false;

    // Check to see if any of the non-muon candidates is inconsistent with a proton 
    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.trk.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // muon should be muon-like...
      if ( idxTrk == idxScdy ) continue; // proton we ID...

      auto const& trk = slc->reco.trk.at(idxTrk);
      if(trk.bestplane == -1) continue; // can't determine chi2 for this track... sigh.

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
				                             slc->vertex.y - trk.start.y,
				                             slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
			                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
			                        	(trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
			                         !isnan(trk.end.y) &&
			                         ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
			                         !isnan(trk.end.z) &&
			                         ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      if ( Atslc >= 10.0 || !trk.pfp.parent_is_primary ) continue; // not a primary so we choose to move on

      // Now check if it's actually proton like and return false if not...
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
      if ( Chi2Proton > 100. || Chi2Muon < 30. ) shouldCutEvt=true;
      else {
        // INFO ON CALO WE WANT TO PRINT OUT...
        const int pdg = trk.truth.p.pdg;
        for ( auto const& point : trk.calo[trk.bestplane].points ) {
          std::cout << pdg << ", " << point.rr << ", " << point.dedx << std::endl;
        }
      }
    }

    if ( shouldCutEvt ) return false;
    // All our primary tracks are proton-like, we're good!
    return true;
  });

const Cut kTertiaryParticleShowers__PRINT([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    int secondaryInd = kPTrackIndNewProton(slc);
    if ( secondaryInd < 0 ) return false;
    unsigned int idxScdy = (unsigned int)secondaryInd;

    // Now check out 'appreciable' showers near the vertex
    for ( auto const& shw : slc->reco.shw ) 
    {
      const float Atslc = std::hypot(slc->vertex.x - shw.start.x,
                                     slc->vertex.y - shw.start.y,
                                     slc->vertex.z - shw.start.z);

      if ( Atslc > 15. || !shw.pfp.parent_is_primary || shw.bestplane_energy < 0.030 ) continue;

      const int pdg = shw.truth.p.pdg;
      std::cout << pdg << ", " << shw.bestplane_energy << ", " << shw.bestplane_dEdx << ", "
                << shw.density << ", " << shw.conversion_gap << ", " << shw.open_angle << std::endl;
    }

    return true;
  });

const Cut kNuMISelection_0Pi__PRINT = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && 
                                        kPTrackNew && kProtonTrack && kNonMuonTracksAreProtons__PRINT );

//const Cut kNuMISelection_0Pi__PRINT = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && 
//                                        kPTrackNew && kProtonTrack && kNonMuonTracksAreProtons && kTertiaryParticleShowers__PRINT );