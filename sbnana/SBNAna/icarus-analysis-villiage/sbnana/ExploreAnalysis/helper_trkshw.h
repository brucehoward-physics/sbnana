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

const Binning kBinsL = Binning::Simple(50,0.,1000.);
const Binning kBinsLE = Binning::Simple(40,0.,4.);

const Binning kBinsCosTh = Binning::Simple(10,-1.,1.);

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
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
    {
        if( slc->reco.pfp.at(i).trackScore < 0.5 ) continue;
        auto const& trk = slc->reco.pfp.at(i).trk;
        if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

        // First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);

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

const Var kRecoMuonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNew(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;
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
        auto const& trk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; // need a different track...
  
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
        TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
        TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      // do we want to make the proton cut even tighter on PID
      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 && trk.len > maxLength ) {
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; // need a different track...
  
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
	      TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
	      TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
	      angle = TMath::Cos(muDir.Angle(pDir));
      }

      // do we want to make the proton cut even tighter on PID
      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && angle >= -0.9 && trk.len > maxLength ) {
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
        auto const& trk = slc->reco.pfp.at(kPTrackIndNewProton(slc)).trk;
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
        auto const& trk = slc->reco.pfp.at(kPTrackIndNewProtonSimple(slc)).trk;
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
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProtonSimple(slc)).trk;
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
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProtonSimple(slc)).trk;
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

// Is selected proton (kPTrackIndNewProton) actually a proton
const Cut kRecoProtonIsTrueProton([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProton(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProton(slc)).trk;

      if ( trk.truth.p.pdg == 2212 ) return true;
      else return false;
    }

    return false;
  });

const Cut kRecoProtonIsTrueProtonSimple([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProtonSimple(slc)).trk;

      if ( trk.truth.p.pdg == 2212 ) return true;
      else return false;
    }

    return false;
  });

const Cut kRecoProtonIsTrueStoppingProtonSimple([](const caf::SRSliceProxy* slc) {
    if ( kPTrackIndNewProtonSimple(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProtonSimple(slc)).trk;

      if ( trk.truth.p.pdg == 2212 && trk.truth.p.daughters.size() == 0 ) return true;
      else return false;
    }

    return false;
  });

const Var kMatchedProtonTrueMomentum([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNewProton(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNewProton(slc)).trk;

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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue;  

      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
        TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
        TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 ) {
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue;

      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
        TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
        TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 ) {
        countP += 1;
	      proton_pfpids.push_back( slc->reco.pfp.at(idxTrk).id );
      }
    }

    for ( unsigned int idxStub = 0; idxStub < slc->reco.stub.size(); ++idxStub ) {
      // Check if stub within 10cm of vertex, contained, and does NOT match one of our already
      // selected proton PFParticles. We'll assume it's a proton if so?

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

      bool matchesProtonPFP = false;
      if ( stub.pfpid >= 0 ) {
        for ( auto const& pfpid : proton_pfpids ) {
          if ( pfpid == stub.pfpid ) matchesProtonPFP = true;
        }
      }

      if ( Atslc < 10.0 && Contained && !matchesProtonPFP ) countP+=1;
    }

    return countP;
  });

// cheated version with some reco
const Var kProtonMult_Reco_Cheated([](const caf::SRSliceProxy* slc) -> int {
    int countP = 0;

    for ( unsigned int idxPfp = 0; idxPfp < slc->reco.pfp.size(); ++idxPfp ) {
      auto const& pfp = slc->reco.pfp.at(idxPfp);
      // I think shower fit is easier? Let's look at the shower fit then:
      if ( pfp.shw.start.x < -5.+std::numeric_limits<double>::epsilon() && 
           pfp.shw.start.x > -5.-std::numeric_limits<double>::epsilon() )
        continue;
      const float Atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                                     slc->vertex.y - pfp.shw.start.y,
                                     slc->vertex.z - pfp.shw.start.z);
      if ( Atslc < 10.0 && pfp.parent_is_primary && pfp.shw.truth.p.pdg == 2212  ) countP+=1;
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
    for ( auto const& pfp : slc.reco.pfp ) {
      // I think shower fit is easier? Let's look at the shower fit then:
      if ( pfp.shw.start.x < -5.+std::numeric_limits<double>::epsilon() && 
           pfp.shw.start.x > -5.-std::numeric_limits<double>::epsilon() )
        continue;
      if ( g4idFound.find( pfp.shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[pfp.shw.truth.p.G4ID]+=1;
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
    for ( auto const& pfp : slc.reco.pfp ) {
      // I think shower fit is easier? Let's look at the shower fit then:
      if ( pfp.shw.start.x < -5.+std::numeric_limits<double>::epsilon() && 
           pfp.shw.start.x > -5.-std::numeric_limits<double>::epsilon() )
        continue;
      if ( g4idFound.find( pfp.shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[pfp.shw.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) protonKEs.push_back( g4idKEs[g4id] );
  }

  return protonKEs;
});

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
    for ( auto const& pfp : slc.reco.pfp ) {
      // I think shower fit is easier? Let's look at the shower fit then:
      if ( pfp.shw.start.x < -5.+std::numeric_limits<double>::epsilon() && 
           pfp.shw.start.x > -5.-std::numeric_limits<double>::epsilon() )
        continue;
      if ( g4idFound.find( pfp.shw.truth.p.G4ID ) != g4idFound.end() ) g4idFound[pfp.shw.truth.p.G4ID]+=1;
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
    for ( auto const& pfp : slc.reco.pfp ) {
      if( pfp.trackScore < 0.5 ) continue;
      auto const& trk = pfp.trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
      // See if it matches our proton G4ID
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idPs[g4id] );
  }

  return ps;
});

// This version does not require it to be track-like, just has a track fit.
const SpillMultiVar kTrueProtonsMomentum_TrackOrShower( [](const caf::SRSpillProxy *sr) {
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
    for ( auto const& pfp : slc.reco.pfp ) {
      auto const& trk = pfp.trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
      // See if it matches our proton G4ID
      if ( g4idFound.find( trk.truth.p.G4ID ) != g4idFound.end() ) g4idFound[trk.truth.p.G4ID]+=1;
    }
  }

  for ( auto const &[g4id, counts] : g4idFound ) {
    if ( counts > 0 ) ps.push_back( g4idPs[g4id] );
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      //if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
      if ( slc->reco.pfp.at(idxTrk).trackScore >= 0.5 && (Chi2Proton > 100. || Chi2Muon < 30.) ) continue;
      if ( slc->reco.pfp.at(idxTrk).trackScore > 0. && slc->reco.pfp.at(idxTrk).trackScore < 0.5 && (Chi2Proton > 100. || Chi2Muon < 40.) ) continue;
      ////////////

      TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
      TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength && TMath::Cos(muDir.Angle(pDir)) >= -0.9 ) {
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      //if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
      //if ( Chi2Proton > 100. || Chi2Muon < 30. ) continue;
      if ( slc->reco.pfp.at(idxTrk).trackScore >= 0.5 && (Chi2Proton > 100. || Chi2Muon < 30.) ) continue;
      if ( slc->reco.pfp.at(idxTrk).trackScore > 0. && slc->reco.pfp.at(idxTrk).trackScore < 0.5 && (Chi2Proton > 100. || Chi2Muon < 40.) ) continue;
      ////////////

      TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
      TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );

      float protonp = -1.;
      if(Contained) protonp = trk.rangeP.p_proton;
      //else {
  	  //  protonp = trk.mcsP.fwdP_proton;
	    //}

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength && TMath::Cos(muDir.Angle(pDir)) >= -0.9 && protonp > 0.4 ) {
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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( primTrks != 2 ) return false;
    if ( maxLength < 0. ) return false;

    TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
    TVector3 piDir( slc->reco.pfp.at(idxScdy).trk.dir.x, slc->reco.pfp.at(idxScdy).trk.dir.y, slc->reco.pfp.at(idxScdy).trk.dir.z );
    if ( TMath::Cos(muDir.Angle(piDir)) < -0.9 ) return false;
    return true;
  }); // kOnlyScdyPionTrack

// No Showers Cut:
const Cut kNoAppreciableShower([](const caf::SRSliceProxy* slc) {
    unsigned int appreciableShowers = 0;
    for ( auto const& pfp : slc->reco.pfp ) 
    {
      if( pfp.trackScore >= 0.5 ) continue;
      auto const& shw = pfp.shw;
      if ( shw.start.x < -5.+std::numeric_limits<double>::epsilon() && 
           shw.start.x > -5.-std::numeric_limits<double>::epsilon() )
        continue;

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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength ) {
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
    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue;

      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      if ( Atslc < 10.0 ) {
        idxScdy = idxTrk;
        break;
      }
    }

    // Muon:
    double pmuon = slc->reco.pfp.at(idxPrim).trk.rangeP.p_muon;
    double Emuon = sqrt((.105658*.105658) + (pmuon*pmuon));
    TVector3 vecpmuon(pmuon*slc->reco.pfp.at(idxPrim).trk.dir.x,pmuon*slc->reco.pfp.at(idxPrim).trk.dir.y,pmuon*slc->reco.pfp.at(idxPrim).trk.dir.z);

    double plmuon = vecpmuon.Dot(NuDirection_NuMI); // longitundinal momentum
    TVector3 ptmuon = vecpmuon.Cross(NuDirection_NuMI); // transverse momentum

    // Pion:
    double ppion = slc->reco.pfp.at(idxScdy).trk.rangeP.p_pion;
    double Epion = sqrt((.139571*.139571) + (ppion*ppion));
    TVector3 vecppion(ppion*slc->reco.pfp.at(idxScdy).trk.dir.x,ppion*slc->reco.pfp.at(idxScdy).trk.dir.y,ppion*slc->reco.pfp.at(idxScdy).trk.dir.z);

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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return false;

    TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
    TVector3 piDir( slc->reco.pfp.at(idxScdy).trk.dir.x, slc->reco.pfp.at(idxScdy).trk.dir.y, slc->reco.pfp.at(idxScdy).trk.dir.z );
    if ( TMath::Cos(muDir.Angle(piDir)) < -0.9 ) return false;
    return true;
  }); // kHasScdyPionTrack

// ***************** *************** *******************
// ********** *************** *************** **********

const Cut kThreePrimaryTracks([](const caf::SRSliceProxy* slc) {
    unsigned int nPrim=0;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary ) nPrim+=1;
    }

    return nPrim==3;
  });

const Cut kNuMIRESSelection = ( kRFiducialNew && kNotClearCosmic && kCutCRLongTrkDirY && kPTrackNew && kHasScdyPionTrack && kNoAppreciableShower && kProtonTrack && kThreePrimaryTracks );

/////// Gray's proposed sample
const Cut kGraysProposedSampleCut([](const caf::SRSliceProxy* slc) {
    if (slc->is_clear_cosmic) return false;
    if (slc->reco.pfp.size() < 3) return false;

    unsigned int nPrimTracks = 1; // think after talking to Jaesung this needs to be 1 - we found this bug while comparing cuts
    bool oneIsProton = false;

    for ( auto const& pfp : slc->reco.pfp ) {
      if( pfp.trackScore < 0.5 ) continue;
      auto const& trk = pfp.trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
    if (slc->reco.pfp.size() < 3) return false;

    unsigned int nPrimTracks = 1; // think after talking to Jaesung this needs to be 1 - we found this bug while comparing cuts
    bool oneIsProton = false;

    int idxMuon = kPTrackIndNew(slc);
    if ( idxMuon < 0 ) return false;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == idxMuon ) continue; // need a different track...
  
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == idxMuon ) continue; // need a different track...

      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

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
    auto const& trk = slc->reco.pfp.at(kProtonTrackIndGraysSample(slc)).trk;
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
    auto const& trk = slc->reco.pfp.at(kProtonTrackIndGraysSample(slc)).trk;

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
          if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 )
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
          if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 )
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
         kProtonTrack400MeV(&slc) ) {
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
          if ( !prim.contained ) continue;
          float p = sqrt(std::pow( prim.startp.x, 2 ) + std::pow( prim.startp.y, 2 ) + std::pow( prim.startp.z, 2 ));
          if ( p > 0.4 )
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
         kProtonTrack400MeV(&slc) ) {
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
         kProtonTrack400MeV(&slc) ) {
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
                                         if ( fabs(thXZ) <= 5. ) continue; // very approximate
                                         thNuMI = dir.Angle(NuDirection_NuMI);
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thNuMI);
                                 });

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
                                         if ( fabs(thXZ) <= 5. ) continue; // very approximate
                                         if ( fabs(thXZ) >= 60. ) continue; // very approximate
                                         thNuMI = dir.Angle(NuDirection_NuMI);
                                       }
                                     }
                                   }

                                   if ( g4id < 0 ) return -999;
                                   return TMath::Cos(thNuMI);
                                 });