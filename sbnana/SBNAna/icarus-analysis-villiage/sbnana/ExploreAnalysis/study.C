// SBNAna
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

//#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

// not in sbn2022a ?
//#include "sbnana/CAFAna/Unfold/UnfoldIterative.h"
//#include "sbnana/CAFAna/Unfold/UnfoldSVD.h"
//#include "sbnana/CAFAna/Unfold/UnfoldTikhonov.h"
//#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Core/Utilities.h"

// don't use -- have mostly updated cuts at moment...
//#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
//#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

// ROOT
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "THStack.h"
#include "TRandom3.h"

// C++
#include <vector>
#include <string>

#include "helper.h"

using namespace ana;

//double pot = 6.0e20;
//double pot = 1.75e20; // Assume POT = 5 months of running with the 0.35e20 POT/month we had during Run 1 -- unfortunately not a very large data set
double pot = 3.0e20; // assume 6e19 per month, then 5 months is 3e20

// NU + COSMIC
// BIGGER DATA SET
const std::string loadstr_ICARUS = "/pnfs/sbn/data/sbn_fd/poms_production/NuMI_Nu_Cosmics_Ovb/mc/reconstructed/icaruscode/v09_37_02_04/flatcaf/*[0,1,2,3,4]*/*/flat*.root";
// FROM JAESUNG
//const std::string loadstr_ICARUS = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_02_04/icarus_numi_nu_cosmics_v09_37_02_04_caf/flat*.root";

// BNB
// Small set
//const std::string loadstr_ICARUS = "/pnfs/sbn/data/sbn_fd/poms_production/BNB_Nu_Cosmics/mc/reconstructed/icaruscode/v09_37_02_04/flatcaf/*0*/*/*flat*.root";

// In-time cosmic:
//const std::string loadstr_ICARUS_intime = "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics_withOverburden2/mc/reconstructed/icaruscode/v09_37_02_07/flatcaf/*0*/*[1,2,3]*/flat*.root";
// FROM JAESUNG
const std::string loadstr_ICARUS_intime = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics_withOverburden2/flatcaf/v09_37_02_07/IcarusProd_2022A_NUMI_in-time_Cosmics_withOverburden_v09_37_02_07_caf/flatcaf*1*.root"; //"[0,1,2,3,4]*.root";

// Slim set
//const std::string loadstr_ICARUS = "/pnfs/sbn/data/sbn_fd/poms_production/NuMI_Nu_Cosmics_Ovb/mc/reconstructed/icaruscode/v09_37_02_04/flatcaf/*0*/*[0,1,2,3,4,5]*/flat*.root";
//const std::string loadstr_ICARUS_intime = "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics_withOverburden2/mc/reconstructed/icaruscode/v09_37_02_07/flatcaf/*0*/*1*/flat*.root";

// VERY few files to test
//const std::string loadstr_ICARUS = "/pnfs/sbn/data/sbn_fd/poms_production/NuMI_Nu_Cosmics_Ovb/mc/reconstructed/icaruscode/v09_37_02_04/flatcaf/*0*/*[0,1,2]*/flat*.root";
//const std::string loadstr_ICARUS_intime = "/pnfs/sbn/data/sbn_fd/poms_production/NUMI_in-time_Cosmics_withOverburden2/mc/reconstructed/icaruscode/v09_37_02_07/flatcaf/*0*/*1*/flat*.root";

const bool fRun = true; //false; //true;

void study ()
{
  if ( fRun ) {
    SpectrumLoader loader( loadstr_ICARUS );
    SpectrumLoader loaderInTime( loadstr_ICARUS_intime );


    // MAKE SPECTRA
    ////////////////////
    Spectrum sAll_NoCut      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNoCut );
    Spectrum sSignal_NoCut   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig );
    Spectrum sOtherNuCC_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll );
    Spectrum sOtherNuNC_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC );
    Spectrum sCosmic_NoCut   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic );
    Spectrum sInTime_NoCut   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNoCut );
    
    Spectrum sAll_Selct      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMISelection );
    Spectrum sSignal_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig && kNuMISelection );
    Spectrum sOtherNuCC_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll && kNuMISelection );
    Spectrum sOtherNuNC_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kNuMISelection );
    Spectrum sCosmic_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kNuMISelection );
    Spectrum sInTime_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNuMISelection );

    // -- Look at the protons in selection
    // kRecoProtonP kRecoProtonIsTrueProton
    Spectrum sAll_P_Selct      ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kNuMISelection );
    Spectrum sSignal_P_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutTrueSig && kNuMISelection );
    Spectrum sOtherNuCC_P_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutNuCCButNotSigAll && kNuMISelection );
    Spectrum sOtherNuNC_P_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutTrueNC && kNuMISelection );
    Spectrum sCosmic_P_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutCosmic && kNuMISelection );
    Spectrum sInTime_P_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loaderInTime, kRecoProtonP, kNoSpillCut, kNuMISelection );

    Spectrum sSignal_TruP_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueProton && kNuMISelection );
    Spectrum sSignal_NotP_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, kCutTrueSig && !kRecoProtonIsTrueProton && kNuMISelection );
    Spectrum sBackgd_P_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonP, kNoSpillCut, !kCutTrueSig && kNuMISelection );

    // kMatchedProtonTrueMomentum
    Spectrum sRecoVsTrue_ProtonP( "True p Momentum [GeV/c]", "Reco p Momentum [GeV/c]",
				  loader,
				  kBinsProtonP, kMatchedProtonTrueMomentum,
				  kBinsProtonP, kRecoProtonP,
				  kNoSpillCut,
				  kCutTrueSig && kRecoProtonIsTrueProton && kNuMISelection );


    // Simple version selecting tracks without Chi2
    Spectrum sAll_PSimple_Selct      ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kNuMISelectionWOProton );
    Spectrum sSignal_PSimple_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutTrueSig && kNuMISelectionWOProton );
    Spectrum sOtherNuCC_PSimple_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutNuCCButNotSigAll && kNuMISelectionWOProton );
    Spectrum sOtherNuNC_PSimple_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutTrueNC && kNuMISelectionWOProton );
    Spectrum sCosmic_PSimple_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutCosmic && kNuMISelectionWOProton );
    Spectrum sInTime_PSimple_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loaderInTime, kRecoProtonPSimple, kNoSpillCut, kNuMISelectionWOProton );

    Spectrum sSignal_TruPSimple_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueProtonSimple && kNuMISelectionWOProton );
    Spectrum sSignal_NotPSimple_Selct( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, kCutTrueSig && !kRecoProtonIsTrueProtonSimple && kNuMISelectionWOProton );
    Spectrum sBackgd_PSimple_Selct   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPSimple, kNoSpillCut, !kCutTrueSig && kNuMISelectionWOProton );

    Spectrum sSignal_TruPSimple_Selct_Chi2P( "Chi2 PID Proton", kBinsPID, loader, kRecoProtonChi2PSimple, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueProtonSimple && kNuMISelectionWOProton );
    Spectrum sSignal_TruPSimple_Selct_Chi2Mu( "Chi2 PID Muon", kBinsPID, loader, kRecoProtonChi2MuSimple, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueProtonSimple && kNuMISelectionWOProton );
    // stopping version:
    Spectrum sSignal_TruStopPSimple_Selct_Chi2P( "Chi2 PID Proton", kBinsPID, loader, kRecoProtonChi2PSimple, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueStoppingProtonSimple && kNuMISelectionWOProton );
    Spectrum sSignal_TruStopPSimple_Selct_Chi2Mu( "Chi2 PID Muon", kBinsPID, loader, kRecoProtonChi2MuSimple, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueStoppingProtonSimple && kNuMISelectionWOProton );


    // and look at multiplicities
    Spectrum sSignal_PMultPreFSI_NoCut      ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PreFSI, kNoSpillCut, kCutTrueSig );
    Spectrum sSignal_PMultPostFSI_NoCut     ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI, kNoSpillCut, kCutTrueSig );
    Spectrum sSignal_PMultPostFSIMinE_NoCut ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI_15MeV, kNoSpillCut, kCutTrueSig );
    // -- truth requiring at least 100MeV
    Spectrum sSignal_PMultPostFSIE100_NoCut ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI_100MeV, kNoSpillCut, kCutTrueSig );

    Spectrum sSignal_PMultPreFSI_Selct      ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PreFSI, kNoSpillCut, kCutTrueSig && kNuMISelection );
    Spectrum sSignal_PMultPostFSI_Selct     ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI, kNoSpillCut, kCutTrueSig && kNuMISelection );
    Spectrum sSignal_PMultPostFSIMinE_Selct ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI_15MeV, kNoSpillCut, kCutTrueSig && kNuMISelection );
    // -- truth requiring at least 100MeV
    Spectrum sSignal_PMultPostFSIE100_Selct ( "Proton Multiplicity", kBinsCounts, loader, kProtonMult_True_PostFSI_100MeV, kNoSpillCut, kCutTrueSig && kNuMISelection );

    // -- Reco versions
    Spectrum sSignal_PMultReco_NoCut        ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco, kNoSpillCut, kCutTrueSig );
    Spectrum sSignal_PMultRecoStubs_NoCut   ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco_wStubs, kNoSpillCut, kCutTrueSig );

    Spectrum sSignal_PMultReco_Selct        ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco, kNoSpillCut, kCutTrueSig && kNuMISelection );
    Spectrum sSignal_PMultRecoStubs_Selct   ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco_wStubs, kNoSpillCut, kCutTrueSig && kNuMISelection );

    // -- cheated version
    Spectrum sSignal_PMultRecoCheated_NoCut ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco_Cheated, kNoSpillCut, kCutTrueSig );
    Spectrum sSignal_PMultRecoCheated_Selct ( "Reco Proton Mult.", kBinsCounts, loader, kProtonMult_Reco_Cheated, kNoSpillCut, kCutTrueSig && kNuMISelection );


    // -- multiplicity and cheated efficiency
    Spectrum sSignal_TrueProtonKE_SpillMV                ( "Proton Kinetic Energy [GeV]", kBinsPZoom, loader, kTrueProtonsKineticEnergy, kNoSpillCut );
    Spectrum sSignal_TrueProtonKE_CheatedReco_SpillMV    ( "Proton Kinetic Energy [GeV]", kBinsPZoom, loader, kTrueProtonsKineticEnergy_CheatedReco, kNoSpillCut );
    Spectrum sSignal_TrueProtonKE_RES_SpillMV            ( "Proton Kinetic Energy [GeV]", kBinsPZoom, loader, kTrueProtonsKineticEnergy_RES, kNoSpillCut );
    Spectrum sSignal_TrueProtonKE_RES_CheatedReco_SpillMV( "Proton Kinetic Energy [GeV]", kBinsPZoom, loader, kTrueProtonsKineticEnergy_RES_CheatedReco, kNoSpillCut );
    // -- momentum version
    Spectrum sSignal_TrueProtonMom_SpillMV                ( "Proton Momentum [GeV/c]", kBinsProtonPZoom, loader, kTrueProtonsMomentum, kNoSpillCut );
    Spectrum sSignal_TrueProtonMom_CheatedReco_SpillMV    ( "Proton Momentum [GeV/c]", kBinsProtonPZoom, loader, kTrueProtonsMomentum_CheatedReco, kNoSpillCut );
    // -- selection efficiency:
    Spectrum sSignal_TrueNuE_SpillMV ( "True #nu Energy [GeV]", kBinsE, loader, kTrueSignalEnu, kNoSpillCut );
    Spectrum sSignal_TrueNuE_Selct_SpillMV ( "True #nu Energy [GeV]", kBinsE, loader, kTrueSignalEnuSelected, kNoSpillCut );
    Spectrum sSignalwProton_TrueMuMom_SpillMV ( "True #mu momentum [GeV/c]", kBinsP, loader, kTrueSignalwProtonMuMom, kNoSpillCut );
    Spectrum sSignalwProton_TrueMuMom_Selct_SpillMV ( "True #mu momentum [GeV/c]", kBinsP, loader, kTrueSignalwProtonMuMomSelected, kNoSpillCut );
    Spectrum sSignalwProton_RecoMuMom_Selct_SpillMV ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kTrueSignalwProtonRecoMuMomSelected, kNoSpillCut );
    Spectrum sAll_RecoMuMom_Selct_SpillMV ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kAllwProtonRecoMuMomSelected, kNoSpillCut );
    Spectrum sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic ( "Reco #mu momentum [GeV/c]", kBinsP, loaderInTime, kAllwProtonRecoMuMomSelected, kNoSpillCut );

    Spectrum sSignal_TrueNuL_SpillMV ( "True Baseline [m]", kBinsL, loader, kTrueSignalLnu, kNoSpillCut );
    Spectrum sSignal_TrueNuLOverE_SpillMV ( "True L/E [m/MeV]", kBinsLE, loader, kTrueSignalLOverEnu, kNoSpillCut );

    Spectrum sNumuCC_TrueNuE_SpillMV ( "True #nu Energy [GeV]", kBinsE, loader, kTrueNumuCCEnu, kNoSpillCut );
    Spectrum sNumuCC_TrueNuL_SpillMV ( "True Baseline [m]", kBinsL, loader, kTrueNumuCCLnu, kNoSpillCut );
    Spectrum sNumuCC_TrueNuLOverE_SpillMV ( "True L/E [m/MeV]", kBinsLE, loader, kTrueNumuCCLOverEnu, kNoSpillCut );


    // ----- Apply same selection as in sAll_RecoMuMom_Selct_SpillMV
    // ----- BASE SELECTION: kNuMI1Mu1PSelection
    // ----- TRUTH GROUPINGS:
    // -------- Signal: k1Mu1P_TrueSigTopology
    // -------- k1Mu1P_SliceIsNumuCC_MuLengthWrong
    // -------- k1Mu1P_SliceIsNumuCC_PWrong
    // -------- k1Mu1P_SliceIsNumuCC_BothWrong
    // -------- k1Mu1P_SliceIsNumuCC_NonFiducial
    // -------- k1Mu1P_SliceIsNumuCC_Other
    // -------- k1Mu1P_SliceIsNuOther
    // -------- k1Mu1P_SliceIsCosmic
    Spectrum sSignalwProton_recoMuMom_BaseSliceSel ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection );
    Spectrum sSignalwProton_recoMuMom_TrueSignal   ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_TrueSigTopology );
    Spectrum sSignalwProton_recoMuMom_NumuMuWrong  ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNumuCC_MuLengthWrong );
    Spectrum sSignalwProton_recoMuMom_NumuPWrong   ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNumuCC_PWrong );
    Spectrum sSignalwProton_recoMuMom_NumuBothWrong( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNumuCC_BothWrong );
    Spectrum sSignalwProton_recoMuMom_NumuNonFidVol( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNumuCC_NonFiducial );
    Spectrum sSignalwProton_recoMuMom_NumuOther    ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNumuCC_Other );
    Spectrum sSignalwProton_recoMuMom_NuOther      ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsNuOther );
    Spectrum sSignalwProton_recoMuMom_Cosmic       ( "Reco #mu momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsCosmic );
    Spectrum sSignalwProton_recoMuMom_Cosmic_InTime( "Reco #mu momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNuMI1Mu1PSelection && k1Mu1P_SliceIsCosmic );

    // -- SpillVars looking at muons from neutrinos somewhat broadly
    Spectrum sNuMIMuons_CosThXZ ( "Cos(ThXZ)", kBinsCosTh, loader, kSingleMatchCosThXZ, kNoSpillCut );
    Spectrum sNuMIMuons_CosThNuMI ( "Cos(ThNuMI)", kBinsCosTh, loader, kSingleMatchCosThNuMI, kNoSpillCut );
    Spectrum sNuMIMuons_CosThNuMINotIsoch ( "Cos(ThNuMI)", kBinsCosTh, loader, kSingleMatchCosThNuMINotIsoch, kNoSpillCut );
    Spectrum sNuMIMuons_CosThNuMINotIsochNotPerp ( "Cos(ThNuMI)", kBinsCosTh, loader, kSingleMatchCosThNuMINotIsochNotPerp, kNoSpillCut );
    

    // Reco w/ stubs vs. PostFSIMinE Truth
    Spectrum sRecoVsTrue_ProtonMult_NoCut( "True p Multiplicity [>15 MeV]", "Reco p Multiplicity [w/ stubs]",
					   loader,
					   kBinsCounts, kProtonMult_True_PostFSI_15MeV,
					   kBinsCounts, kProtonMult_Reco_wStubs,
					   kNoSpillCut,
					   kCutTrueSig );

    Spectrum sRecoVsTrue_ProtonMult_Selct( "True p Multiplicity [>15 MeV]", "Reco p Multiplicity [w/ stubs]",
                                           loader,
                                           kBinsCounts, kProtonMult_True_PostFSI_15MeV,
                                           kBinsCounts, kProtonMult_Reco_wStubs,
                                           kNoSpillCut,
                                           kCutTrueSig && kNuMISelection );


    // Interaction Modes
    Spectrum sSignalQEL_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL && kNuMISelection );
    Spectrum sSignalMEC_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC && kNuMISelection );
    Spectrum sSignalRES_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES && kNuMISelection );
    Spectrum sSignalDIS_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS && kNuMISelection );
    Spectrum sSignalCOH_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH && kNuMISelection );
    
    // Contained versions
    Spectrum sSignal_Cont_NoCut   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigContained );
    Spectrum sOtherNuCC_Cont_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigContained );
    
    Spectrum sAll_Cont_Selct      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMISelectionContained );
    Spectrum sSignal_Cont_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigContained && kNuMISelectionContained );
    Spectrum sOtherNuCC_Cont_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigContained && kNuMISelectionContained );
    Spectrum sOtherNuNC_Cont_Selct( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kNuMISelectionContained );
    Spectrum sCosmic_Cont_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kNuMISelectionContained );
    Spectrum sInTime_Cont_Selct   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNuMISelectionContained );
    
    // Gray's sample
    Spectrum sAll_SelctG      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignal_SelctG   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sOtherNuCC_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sOtherNuNC_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sCosmic_SelctG   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sInTime_SelctG   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack );

    // Gray's sample broken into interaction mode
    Spectrum sSignalQEL_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignalMEC_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignalRES_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignalDIS_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignalCOH_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignalELS_SelctG( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigELS && kGraysProposedSampleCutWithMuonTrack );

    // GRAY's selection with proton vars
    // kProtonTrackIndGraysSample,  kRecoProtonPGraysSample,  kRecoProtonIsTrueProtonGraysSample
    Spectrum sAll_P_SelctG      ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignal_P_SelctG   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutTrueSig && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sOtherNuCC_P_SelctG( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutNuCCButNotSigAll && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sOtherNuNC_P_SelctG( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutTrueNC && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sCosmic_P_SelctG   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutCosmic && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sInTime_P_SelctG   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loaderInTime, kRecoProtonPGraysSample, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack );

    Spectrum sSignal_TruP_SelctG( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutTrueSig && kRecoProtonIsTrueProtonGraysSample && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sSignal_NotP_SelctG( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, kCutTrueSig && !kRecoProtonIsTrueProtonGraysSample && kGraysProposedSampleCutWithMuonTrack );
    Spectrum sBackgd_P_SelctG   ( "Reco p Momentum [GeV/c]", kBinsProtonP, loader, kRecoProtonPGraysSample, kNoSpillCut, !kCutTrueSig && kGraysProposedSampleCutWithMuonTrack );



    // Gray's sample with full sel?
    Spectrum sAll_SelTot      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignal_SelTot   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sOtherNuCC_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sOtherNuNC_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sCosmic_SelTot   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sInTime_SelTot   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalQEL_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalMEC_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalRES_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalDIS_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalCOH_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );
    Spectrum sSignalELS_SelTot( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigELS && kGraysProposedSampleCutWithMuonTrack && kNuMISelection );

    // LOOK AT A COHERENT ENHANCING SELECTION?
    // We have the ALL, NC, NonSignal CC, and Cosmic, just need the interaction breakdowns
    Spectrum sSignalQEL_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL );
    Spectrum sSignalMEC_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC );
    Spectrum sSignalRES_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES );
    Spectrum sSignalDIS_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS );
    Spectrum sSignalCOH_NoCut( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH );
    // Selection kNuMICOHSelection
    Spectrum sAll_SelCoh      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMICOHSelection );
    Spectrum sSignal_SelCoh   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig && kNuMICOHSelection );
    Spectrum sOtherNuCC_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll && kNuMICOHSelection );
    Spectrum sOtherNuNC_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kNuMICOHSelection );
    Spectrum sCosmic_SelCoh   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kNuMICOHSelection );
    Spectrum sInTime_SelCoh   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNuMICOHSelection );
    Spectrum sSignalQEL_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL && kNuMICOHSelection );
    Spectrum sSignalMEC_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC && kNuMICOHSelection );
    Spectrum sSignalRES_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES && kNuMICOHSelection );
    Spectrum sSignalDIS_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS && kNuMICOHSelection );
    Spectrum sSignalCOH_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH && kNuMICOHSelection );
    Spectrum sSignalELS_SelCoh( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigELS && kNuMICOHSelection );

    // kTCoherent
    Spectrum sT_All_SelCoh      ( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kNuMICOHSelection );
    Spectrum sT_Signal_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSig && kNuMICOHSelection );
    Spectrum sT_OtherNuCC_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutNuCCButNotSigAll && kNuMICOHSelection );
    Spectrum sT_OtherNuNC_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueNC && kNuMICOHSelection );
    Spectrum sT_Cosmic_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutCosmic && kNuMICOHSelection );
    Spectrum sT_InTime_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsT, loaderInTime, kTCoherent, kNoSpillCut, kNuMICOHSelection );
    Spectrum sT_SignalQEL_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigQEL && kNuMICOHSelection );
    Spectrum sT_SignalMEC_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigMEC && kNuMICOHSelection );
    Spectrum sT_SignalRES_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigRES && kNuMICOHSelection );
    Spectrum sT_SignalDIS_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigDIS && kNuMICOHSelection );
    Spectrum sT_SignalCOH_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigCOH && kNuMICOHSelection );
    Spectrum sT_SignalELS_SelCoh( "Reco |t| [GeV/c]^2", kBinsT, loader, kTCoherent, kNoSpillCut, kCutTrueSigELS && kNuMICOHSelection );

    Spectrum sTzoom_All_SelCoh      ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kNuMICOHSelection );
    Spectrum sTzoom_Signal_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSig && kNuMICOHSelection );
    Spectrum sTzoom_OtherNuCC_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutNuCCButNotSigAll && kNuMICOHSelection );
    Spectrum sTzoom_OtherNuNC_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueNC && kNuMICOHSelection );
    Spectrum sTzoom_Cosmic_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutCosmic && kNuMICOHSelection );
    Spectrum sTzoom_InTime_SelCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loaderInTime, kTCoherent, kNoSpillCut, kNuMICOHSelection );
    Spectrum sTzoom_SignalQEL_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigQEL && kNuMICOHSelection );
    Spectrum sTzoom_SignalMEC_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigMEC && kNuMICOHSelection );
    Spectrum sTzoom_SignalRES_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigRES && kNuMICOHSelection );
    Spectrum sTzoom_SignalDIS_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigDIS && kNuMICOHSelection );
    Spectrum sTzoom_SignalCOH_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigCOH && kNuMICOHSelection );
    Spectrum sTzoom_SignalELS_SelCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigELS && kNuMICOHSelection );

    // with additional cheated cut on proton (40 MeV)
    Spectrum sTzoom_All_ChtCoh      ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_Signal_ChtCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSig && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_OtherNuCC_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutNuCCButNotSigAll && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_OtherNuNC_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueNC && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_Cosmic_ChtCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutCosmic && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_InTime_ChtCoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loaderInTime, kTCoherent, kNoSpillCut, kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalQEL_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigQEL && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalMEC_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigMEC && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalRES_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigRES && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalDIS_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigDIS && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalCOH_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigCOH && kNuMICOHSelection_CheatedNoP );
    Spectrum sTzoom_SignalELS_ChtCoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigELS && kNuMICOHSelection_CheatedNoP );

    // with additional cheated cut on proton (15 MeV)
    Spectrum sTzoom_All_CLECoh      ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_Signal_CLECoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSig && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_OtherNuCC_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutNuCCButNotSigAll && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_OtherNuNC_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueNC && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_Cosmic_CLECoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutCosmic && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_InTime_CLECoh   ( "Reco |t| [GeV/c]^2", kBinsTzoom, loaderInTime, kTCoherent, kNoSpillCut, kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalQEL_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigQEL && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalMEC_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigMEC && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalRES_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigRES && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalDIS_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigDIS && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalCOH_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigCOH && kNuMICOHSelection_CheatedNoPLowE );
    Spectrum sTzoom_SignalELS_CLECoh( "Reco |t| [GeV/c]^2", kBinsTzoom, loader, kTCoherent, kNoSpillCut, kCutTrueSigELS && kNuMICOHSelection_CheatedNoPLowE );



    // RES Selection
    Spectrum sAll_SelRes      ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kNuMIRESSelection );
    Spectrum sSignal_SelRes   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSig && kNuMIRESSelection );
    Spectrum sOtherNuCC_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutNuCCButNotSigAll && kNuMIRESSelection );
    Spectrum sOtherNuNC_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueNC && kNuMIRESSelection );
    Spectrum sCosmic_SelRes   ( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutCosmic && kNuMIRESSelection );
    Spectrum sInTime_SelRes   ( "Reco Momentum [GeV/c]", kBinsP, loaderInTime, kRecoMuonPNew, kNoSpillCut, kNuMIRESSelection );
    Spectrum sSignalQEL_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigQEL && kNuMIRESSelection );
    Spectrum sSignalMEC_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigMEC && kNuMIRESSelection );
    Spectrum sSignalRES_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigRES && kNuMIRESSelection );
    Spectrum sSignalDIS_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigDIS && kNuMIRESSelection );
    Spectrum sSignalCOH_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigCOH && kNuMIRESSelection );
    Spectrum sSignalELS_SelRes( "Reco Momentum [GeV/c]", kBinsP, loader, kRecoMuonPNew, kNoSpillCut, kCutTrueSigELS && kNuMIRESSelection );


    // GO!
    ////////////////////
    loader.Go();
    loaderInTime.Go();

    // Save spectra
    ////////////////////
    TFile *fSpec = new TFile("spectra.root","RECREATE");

    sAll_NoCut.SaveTo(       fSpec->mkdir("sAll_NoCut") );
    sSignal_NoCut.SaveTo(    fSpec->mkdir("sSignal_NoCut") );
    sOtherNuCC_NoCut.SaveTo( fSpec->mkdir("sOtherNuCC_NoCut") );
    sOtherNuNC_NoCut.SaveTo( fSpec->mkdir("sOtherNuNC_NoCut") );
    sCosmic_NoCut.SaveTo(    fSpec->mkdir("sCosmic_NoCut") );
    sInTime_NoCut.SaveTo(    fSpec->mkdir("sInTime_NoCut") );

    sSignal_Cont_NoCut.SaveTo(    fSpec->mkdir("sSignal_Cont_NoCut") );
    sOtherNuCC_Cont_NoCut.SaveTo( fSpec->mkdir("sOtherNuCC_Cont_NoCut") );

    sAll_Selct.SaveTo(       fSpec->mkdir("sAll_Selct") );
    sSignal_Selct.SaveTo(    fSpec->mkdir("sSignal_Selct") );
    sOtherNuCC_Selct.SaveTo( fSpec->mkdir("sOtherNuCC_Selct") );
    sOtherNuNC_Selct.SaveTo( fSpec->mkdir("sOtherNuNC_Selct") );
    sCosmic_Selct.SaveTo(    fSpec->mkdir("sCosmic_Selct") );
    sInTime_Selct.SaveTo(    fSpec->mkdir("sInTime_Selct") );

    sAll_P_Selct.SaveTo(       fSpec->mkdir("sAll_P_Selct") );
    sSignal_P_Selct.SaveTo(    fSpec->mkdir("sSignal_P_Selct") );
    sOtherNuCC_P_Selct.SaveTo( fSpec->mkdir("sOtherNuCC_P_Selct") );
    sOtherNuNC_P_Selct.SaveTo( fSpec->mkdir("sOtherNuNC_P_Selct") );
    sCosmic_P_Selct.SaveTo(    fSpec->mkdir("sCosmic_P_Selct") );
    sInTime_P_Selct.SaveTo(    fSpec->mkdir("sInTime_P_Selct") );
    sSignal_TruP_Selct.SaveTo( fSpec->mkdir("sSignal_TruP_Selct") );
    sSignal_NotP_Selct.SaveTo( fSpec->mkdir("sSignal_NotP_Selct") );
    sBackgd_P_Selct.SaveTo(    fSpec->mkdir("sBackgd_P_Selct") );

    sRecoVsTrue_ProtonP.SaveTo( fSpec->mkdir("sRecoVsTrue_ProtonP") );


    sAll_PSimple_Selct.SaveTo( fSpec->mkdir("sAll_PSimple_Selct") );
    sSignal_PSimple_Selct.SaveTo( fSpec->mkdir("sSignal_PSimple_Selct") );
    sOtherNuCC_PSimple_Selct.SaveTo( fSpec->mkdir("sOtherNuCC_PSimple_Selct") );
    sOtherNuNC_PSimple_Selct.SaveTo( fSpec->mkdir("sOtherNuNC_PSimple_Selct") );
    sCosmic_PSimple_Selct.SaveTo( fSpec->mkdir("sCosmic_PSimple_Selct") );
    sInTime_PSimple_Selct.SaveTo( fSpec->mkdir("sInTime_PSimple_Selct") );

    sSignal_TruPSimple_Selct.SaveTo( fSpec->mkdir("sSignal_TruPSimple_Selct") );
    sSignal_NotPSimple_Selct.SaveTo( fSpec->mkdir("sSignal_NotPSimple_Selct") );
    sBackgd_PSimple_Selct.SaveTo( fSpec->mkdir("sBackgd_PSimple_Selct") );

    sSignal_TruPSimple_Selct_Chi2P.SaveTo( fSpec->mkdir("sSignal_TruPSimple_Selct_Chi2P") );
    sSignal_TruPSimple_Selct_Chi2Mu.SaveTo( fSpec->mkdir("sSignal_TruPSimple_Selct_Chi2Mu") );
    sSignal_TruStopPSimple_Selct_Chi2P.SaveTo( fSpec->mkdir("sSignal_TruStopPSimple_Selct_Chi2P") );
    sSignal_TruStopPSimple_Selct_Chi2Mu.SaveTo( fSpec->mkdir("sSignal_TruStopPSimple_Selct_Chi2Mu") );



    sSignal_PMultPreFSI_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultPreFSI_NoCut") );
    sSignal_PMultPostFSI_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultPostFSI_NoCut") );
    sSignal_PMultPostFSIMinE_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultPostFSIMinE_NoCut") );
    sSignal_PMultReco_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultReco_NoCut") );
    sSignal_PMultRecoStubs_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultRecoStubs_NoCut") );

    sSignal_PMultPreFSI_Selct.SaveTo( fSpec->mkdir("sSignal_PMultPreFSI_Selct") );
    sSignal_PMultPostFSI_Selct.SaveTo( fSpec->mkdir("sSignal_PMultPostFSI_Selct") );
    sSignal_PMultPostFSIMinE_Selct.SaveTo( fSpec->mkdir("sSignal_PMultPostFSIMinE_Selct") );
    sSignal_PMultReco_Selct.SaveTo( fSpec->mkdir("sSignal_PMultReco_Selct") );
    sSignal_PMultRecoStubs_Selct.SaveTo( fSpec->mkdir("sSignal_PMultRecoStubs_Selct") );


    sSignal_PMultPostFSIE100_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultPostFSIE100_NoCut") );
    sSignal_PMultPostFSIE100_Selct.SaveTo( fSpec->mkdir("sSignal_PMultPostFSIE100_Selct") );
    sSignal_PMultRecoCheated_NoCut.SaveTo( fSpec->mkdir("sSignal_PMultRecoCheated_NoCut") );
    sSignal_PMultRecoCheated_Selct.SaveTo( fSpec->mkdir("sSignal_PMultRecoCheated_Selct") );

    sSignal_TrueProtonKE_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonKE_SpillMV") );
    sSignal_TrueProtonKE_RES_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonKE_RES_SpillMV") );
    sSignal_TrueProtonKE_CheatedReco_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonKE_CheatedReco_SpillMV") );
    sSignal_TrueProtonKE_RES_CheatedReco_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonKE_RES_CheatedReco_SpillMV") );

    sSignal_TrueProtonMom_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonMom_SpillMV") );
    sSignal_TrueProtonMom_CheatedReco_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueProtonMom_CheatedReco_SpillMV") );

    sSignal_TrueNuE_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueNuE_SpillMV") );
    sSignal_TrueNuE_Selct_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueNuE_Selct_SpillMV") );
    sSignalwProton_TrueMuMom_SpillMV.SaveTo( fSpec->mkdir("sSignalwProton_TrueMuMom_SpillMV") );
    sSignalwProton_TrueMuMom_Selct_SpillMV.SaveTo( fSpec->mkdir("sSignalwProton_TrueMuMom_Selct_SpillMV") );
    sSignalwProton_RecoMuMom_Selct_SpillMV.SaveTo( fSpec->mkdir("sSignalwProton_RecoMuMom_Selct_SpillMV") );
    sAll_RecoMuMom_Selct_SpillMV.SaveTo( fSpec->mkdir("sAll_RecoMuMom_Selct_SpillMV") );
    sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic.SaveTo( fSpec->mkdir("sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic") );

    sSignalwProton_recoMuMom_BaseSliceSel.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_BaseSliceSel") );
    sSignalwProton_recoMuMom_TrueSignal.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_TrueSignal") );
    sSignalwProton_recoMuMom_NumuMuWrong.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NumuMuWrong") );
    sSignalwProton_recoMuMom_NumuPWrong.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NumuPWrong") );
    sSignalwProton_recoMuMom_NumuBothWrong.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NumuBothWrong") );
    sSignalwProton_recoMuMom_NumuNonFidVol.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NumuNonFidVol") );
    sSignalwProton_recoMuMom_NumuOther.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NumuOther") );
    sSignalwProton_recoMuMom_NuOther.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_NuOther") );
    sSignalwProton_recoMuMom_Cosmic.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_Cosmic") );
    sSignalwProton_recoMuMom_Cosmic_InTime.SaveTo( fSpec->mkdir("sSignalwProton_recoMuMom_Cosmic_InTime") );

    sSignal_TrueNuL_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueNuL_SpillMV") );
    sSignal_TrueNuLOverE_SpillMV.SaveTo( fSpec->mkdir("sSignal_TrueNuLOverE_SpillMV") );

    sNumuCC_TrueNuE_SpillMV.SaveTo( fSpec->mkdir("sNumuCC_TrueNuE_SpillMV") );    
    sNumuCC_TrueNuL_SpillMV.SaveTo( fSpec->mkdir("sNumuCC_TrueNuL_SpillMV") );
    sNumuCC_TrueNuLOverE_SpillMV.SaveTo( fSpec->mkdir("sNumuCC_TrueNuLOverE_SpillMV") );


    sNuMIMuons_CosThXZ.SaveTo( fSpec->mkdir("sNuMIMuons_CosThXZ") );
    sNuMIMuons_CosThNuMI.SaveTo( fSpec->mkdir("sNuMIMuons_CosThNuMI") );
    sNuMIMuons_CosThNuMINotIsoch.SaveTo( fSpec->mkdir("sNuMIMuons_CosThNuMINotIsoch") );
    sNuMIMuons_CosThNuMINotIsochNotPerp.SaveTo( fSpec->mkdir("sNuMIMuons_CosThNuMINotIsochNotPerp") );



    sRecoVsTrue_ProtonMult_NoCut.SaveTo( fSpec->mkdir("sRecoVsTrue_ProtonMult_NoCut") );
    sRecoVsTrue_ProtonMult_Selct.SaveTo( fSpec->mkdir("sRecoVsTrue_ProtonMult_Selct") );

    sSignalQEL_Selct.SaveTo( fSpec->mkdir("sSignalQEL_Selct") );
    sSignalMEC_Selct.SaveTo( fSpec->mkdir("sSignalMEC_Selct") );
    sSignalRES_Selct.SaveTo( fSpec->mkdir("sSignalRES_Selct") );
    sSignalDIS_Selct.SaveTo( fSpec->mkdir("sSignalDIS_Selct") );
    sSignalCOH_Selct.SaveTo( fSpec->mkdir("sSignalCOH_Selct") );

    sAll_Cont_Selct.SaveTo(       fSpec->mkdir("sAll_Cont_Selct") );
    sSignal_Cont_Selct.SaveTo(    fSpec->mkdir("sSignal_Cont_Selct") );
    sOtherNuCC_Cont_Selct.SaveTo( fSpec->mkdir("sOtherNuCC_Cont_Selct") );
    sOtherNuNC_Cont_Selct.SaveTo( fSpec->mkdir("sOtherNuNC_Cont_Selct") );
    sCosmic_Cont_Selct.SaveTo(    fSpec->mkdir("sCosmic_Cont_Selct") );
    sInTime_Cont_Selct.SaveTo(    fSpec->mkdir("sInTime_Cont_Selct") );

    sAll_SelctG.SaveTo(       fSpec->mkdir("sAll_SelctG") );
    sSignal_SelctG.SaveTo(    fSpec->mkdir("sSignal_SelctG") );
    sOtherNuCC_SelctG.SaveTo( fSpec->mkdir("sOtherNuCC_SelctG") );
    sOtherNuNC_SelctG.SaveTo( fSpec->mkdir("sOtherNuNC_SelctG") );
    sCosmic_SelctG.SaveTo(    fSpec->mkdir("sCosmic_SelctG") );
    sInTime_SelctG.SaveTo(    fSpec->mkdir("sInTime_SelctG") );

    sSignalQEL_SelctG.SaveTo( fSpec->mkdir("sSignalQEL_SelctG") );
    sSignalMEC_SelctG.SaveTo( fSpec->mkdir("sSignalMEC_SelctG") );
    sSignalRES_SelctG.SaveTo( fSpec->mkdir("sSignalRES_SelctG") );
    sSignalDIS_SelctG.SaveTo( fSpec->mkdir("sSignalDIS_SelctG") );
    sSignalCOH_SelctG.SaveTo( fSpec->mkdir("sSignalCOH_SelctG") );
    sSignalELS_SelctG.SaveTo( fSpec->mkdir("sSignalELS_SelctG") );


    sAll_P_SelctG.SaveTo(       fSpec->mkdir("sAll_P_SelctG") );
    sSignal_P_SelctG.SaveTo(    fSpec->mkdir("sSignal_P_SelctG") );
    sOtherNuCC_P_SelctG.SaveTo( fSpec->mkdir("sOtherNuCC_P_SelctG") );
    sOtherNuNC_P_SelctG.SaveTo( fSpec->mkdir("sOtherNuNC_P_SelctG") );
    sCosmic_P_SelctG.SaveTo(    fSpec->mkdir("sCosmic_P_SelctG") );
    sInTime_P_SelctG.SaveTo(    fSpec->mkdir("sInTime_P_SelctG") );
    sSignal_TruP_SelctG.SaveTo( fSpec->mkdir("sSignal_TruP_SelctG") );
    sSignal_NotP_SelctG.SaveTo( fSpec->mkdir("sSignal_NotP_SelctG") );
    sBackgd_P_SelctG.SaveTo(    fSpec->mkdir("sBackgd_P_SelctG") );


    sAll_SelTot.SaveTo(       fSpec->mkdir("sAll_SelTot") );
    sSignal_SelTot.SaveTo(    fSpec->mkdir("sSignal_SelTot") );
    sOtherNuCC_SelTot.SaveTo( fSpec->mkdir("sOtherNuCC_SelTot") );
    sOtherNuNC_SelTot.SaveTo( fSpec->mkdir("sOtherNuNC_SelTot") );
    sCosmic_SelTot.SaveTo(    fSpec->mkdir("sCosmic_SelTot") );
    sInTime_SelTot.SaveTo(    fSpec->mkdir("sInTime_SelTot") );
    sSignalQEL_SelTot.SaveTo( fSpec->mkdir("sSignalQEL_SelTot") );
    sSignalMEC_SelTot.SaveTo( fSpec->mkdir("sSignalMEC_SelTot") );
    sSignalRES_SelTot.SaveTo( fSpec->mkdir("sSignalRES_SelTot") );
    sSignalDIS_SelTot.SaveTo( fSpec->mkdir("sSignalDIS_SelTot") );
    sSignalCOH_SelTot.SaveTo( fSpec->mkdir("sSignalCOH_SelTot") );
    sSignalELS_SelTot.SaveTo( fSpec->mkdir("sSignalELS_SelTot") );

    sSignalQEL_NoCut.SaveTo( fSpec->mkdir("sSignalQEL_NoCut") );
    sSignalMEC_NoCut.SaveTo( fSpec->mkdir("sSignalMEC_NoCut") );
    sSignalRES_NoCut.SaveTo( fSpec->mkdir("sSignalRES_NoCut") );
    sSignalDIS_NoCut.SaveTo( fSpec->mkdir("sSignalDIS_NoCut") );
    sSignalCOH_NoCut.SaveTo( fSpec->mkdir("sSignalCOH_NoCut") );
    sAll_SelCoh.SaveTo(       fSpec->mkdir("sAll_SelCoh") );
    sSignal_SelCoh.SaveTo(    fSpec->mkdir("sSignal_SelCoh") );
    sOtherNuCC_SelCoh.SaveTo( fSpec->mkdir("sOtherNuCC_SelCoh") );
    sOtherNuNC_SelCoh.SaveTo( fSpec->mkdir("sOtherNuNC_SelCoh") );
    sCosmic_SelCoh.SaveTo(    fSpec->mkdir("sCosmic_SelCoh") );
    sInTime_SelCoh.SaveTo(    fSpec->mkdir("sInTime_SelCoh") );
    sSignalQEL_SelCoh.SaveTo( fSpec->mkdir("sSignalQEL_SelCoh") );
    sSignalMEC_SelCoh.SaveTo( fSpec->mkdir("sSignalMEC_SelCoh") );
    sSignalRES_SelCoh.SaveTo( fSpec->mkdir("sSignalRES_SelCoh") );
    sSignalDIS_SelCoh.SaveTo( fSpec->mkdir("sSignalDIS_SelCoh") );
    sSignalCOH_SelCoh.SaveTo( fSpec->mkdir("sSignalCOH_SelCoh") );

    sT_All_SelCoh.SaveTo(       fSpec->mkdir("sT_All_SelCoh") );
    sT_Signal_SelCoh.SaveTo(    fSpec->mkdir("sT_Signal_SelCoh") );
    sT_OtherNuCC_SelCoh.SaveTo( fSpec->mkdir("sT_OtherNuCC_SelCoh") );
    sT_OtherNuNC_SelCoh.SaveTo( fSpec->mkdir("sT_OtherNuNC_SelCoh") );
    sT_Cosmic_SelCoh.SaveTo(    fSpec->mkdir("sT_Cosmic_SelCoh") );
    sT_InTime_SelCoh.SaveTo(    fSpec->mkdir("sT_InTime_SelCoh") );
    sT_SignalQEL_SelCoh.SaveTo( fSpec->mkdir("sT_SignalQEL_SelCoh") );
    sT_SignalMEC_SelCoh.SaveTo( fSpec->mkdir("sT_SignalMEC_SelCoh") );
    sT_SignalRES_SelCoh.SaveTo( fSpec->mkdir("sT_SignalRES_SelCoh") );
    sT_SignalDIS_SelCoh.SaveTo( fSpec->mkdir("sT_SignalDIS_SelCoh") );
    sT_SignalCOH_SelCoh.SaveTo( fSpec->mkdir("sT_SignalCOH_SelCoh") );

    sTzoom_All_SelCoh.SaveTo(       fSpec->mkdir("sTzoom_All_SelCoh") );
    sTzoom_Signal_SelCoh.SaveTo(    fSpec->mkdir("sTzoom_Signal_SelCoh") );
    sTzoom_OtherNuCC_SelCoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuCC_SelCoh") );
    sTzoom_OtherNuNC_SelCoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuNC_SelCoh") );
    sTzoom_Cosmic_SelCoh.SaveTo(    fSpec->mkdir("sTzoom_Cosmic_SelCoh") );
    sTzoom_InTime_SelCoh.SaveTo(    fSpec->mkdir("sTzoom_InTime_SelCoh") );
    sTzoom_SignalQEL_SelCoh.SaveTo( fSpec->mkdir("sTzoom_SignalQEL_SelCoh") );
    sTzoom_SignalMEC_SelCoh.SaveTo( fSpec->mkdir("sTzoom_SignalMEC_SelCoh") );
    sTzoom_SignalRES_SelCoh.SaveTo( fSpec->mkdir("sTzoom_SignalRES_SelCoh") );
    sTzoom_SignalDIS_SelCoh.SaveTo( fSpec->mkdir("sTzoom_SignalDIS_SelCoh") );
    sTzoom_SignalCOH_SelCoh.SaveTo( fSpec->mkdir("sTzoom_SignalCOH_SelCoh") );

    sTzoom_All_ChtCoh.SaveTo(       fSpec->mkdir("sTzoom_All_ChtCoh") );
    sTzoom_Signal_ChtCoh.SaveTo(    fSpec->mkdir("sTzoom_Signal_ChtCoh") );
    sTzoom_OtherNuCC_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuCC_ChtCoh") );
    sTzoom_OtherNuNC_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuNC_ChtCoh") );
    sTzoom_Cosmic_ChtCoh.SaveTo(    fSpec->mkdir("sTzoom_Cosmic_ChtCoh") );
    sTzoom_InTime_ChtCoh.SaveTo(    fSpec->mkdir("sTzoom_InTime_ChtCoh") );
    sTzoom_SignalQEL_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_SignalQEL_ChtCoh") );
    sTzoom_SignalMEC_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_SignalMEC_ChtCoh") );
    sTzoom_SignalRES_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_SignalRES_ChtCoh") );
    sTzoom_SignalDIS_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_SignalDIS_ChtCoh") );
    sTzoom_SignalCOH_ChtCoh.SaveTo( fSpec->mkdir("sTzoom_SignalCOH_ChtCoh") );
    //
    sTzoom_All_CLECoh.SaveTo(       fSpec->mkdir("sTzoom_All_CLECoh") );
    sTzoom_Signal_CLECoh.SaveTo(    fSpec->mkdir("sTzoom_Signal_CLECoh") );
    sTzoom_OtherNuCC_CLECoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuCC_CLECoh") );
    sTzoom_OtherNuNC_CLECoh.SaveTo( fSpec->mkdir("sTzoom_OtherNuNC_CLECoh") );
    sTzoom_Cosmic_CLECoh.SaveTo(    fSpec->mkdir("sTzoom_Cosmic_CLECoh") );
    sTzoom_InTime_CLECoh.SaveTo(    fSpec->mkdir("sTzoom_InTime_CLECoh") );
    sTzoom_SignalQEL_CLECoh.SaveTo( fSpec->mkdir("sTzoom_SignalQEL_CLECoh") );
    sTzoom_SignalMEC_CLECoh.SaveTo( fSpec->mkdir("sTzoom_SignalMEC_CLECoh") );
    sTzoom_SignalRES_CLECoh.SaveTo( fSpec->mkdir("sTzoom_SignalRES_CLECoh") );
    sTzoom_SignalDIS_CLECoh.SaveTo( fSpec->mkdir("sTzoom_SignalDIS_CLECoh") );
    sTzoom_SignalCOH_CLECoh.SaveTo( fSpec->mkdir("sTzoom_SignalCOH_CLECoh") );

    sAll_SelRes.SaveTo(       fSpec->mkdir("sAll_SelRes") );
    sSignal_SelRes.SaveTo(    fSpec->mkdir("sSignal_SelRes") );
    sOtherNuCC_SelRes.SaveTo( fSpec->mkdir("sOtherNuCC_SelRes") );
    sOtherNuNC_SelRes.SaveTo( fSpec->mkdir("sOtherNuNC_SelRes") );
    sCosmic_SelRes.SaveTo(    fSpec->mkdir("sCosmic_SelRes") );
    sInTime_SelRes.SaveTo(    fSpec->mkdir("sInTime_SelRes") );
    sSignalQEL_SelRes.SaveTo( fSpec->mkdir("sSignalQEL_SelRes") );
    sSignalMEC_SelRes.SaveTo( fSpec->mkdir("sSignalMEC_SelRes") );
    sSignalRES_SelRes.SaveTo( fSpec->mkdir("sSignalRES_SelRes") );
    sSignalDIS_SelRes.SaveTo( fSpec->mkdir("sSignalDIS_SelRes") );
    sSignalCOH_SelRes.SaveTo( fSpec->mkdir("sSignalCOH_SelRes") );

    fSpec->Close();
  }

  // HISTOGRAMS
  ////////////////////
  std::string fLoad = "spectra.root";

  Spectrum *sAll_NoCut = LoadFromFile<Spectrum>(fLoad,"sAll_NoCut").release();
  Spectrum *sSignal_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_NoCut").release();
  Spectrum *sOtherNuCC_NoCut = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_NoCut").release();
  Spectrum *sOtherNuNC_NoCut = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_NoCut").release();
  Spectrum *sCosmic_NoCut = LoadFromFile<Spectrum>(fLoad,"sCosmic_NoCut").release();
  Spectrum *sInTime_NoCut = LoadFromFile<Spectrum>(fLoad,"sInTime_NoCut").release();

  Spectrum *sSignalQEL_Selct = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_Selct").release();
  Spectrum *sSignalMEC_Selct = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_Selct").release();
  Spectrum *sSignalRES_Selct = LoadFromFile<Spectrum>(fLoad,"sSignalRES_Selct").release();
  Spectrum *sSignalDIS_Selct = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_Selct").release();
  Spectrum *sSignalCOH_Selct = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_Selct").release();

  Spectrum *sSignal_Cont_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_Cont_NoCut").release();
  Spectrum *sOtherNuCC_Cont_NoCut = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_Cont_NoCut").release();

  Spectrum *sAll_Selct = LoadFromFile<Spectrum>(fLoad,"sAll_Selct").release();
  Spectrum *sSignal_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_Selct").release();
  Spectrum *sOtherNuCC_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_Selct").release();
  Spectrum *sOtherNuNC_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_Selct").release();
  Spectrum *sCosmic_Selct = LoadFromFile<Spectrum>(fLoad,"sCosmic_Selct").release();
  Spectrum *sInTime_Selct = LoadFromFile<Spectrum>(fLoad,"sInTime_Selct").release();

  Spectrum *sAll_P_Selct = LoadFromFile<Spectrum>(fLoad,"sAll_P_Selct").release();
  Spectrum *sSignal_P_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_P_Selct").release();
  Spectrum *sOtherNuCC_P_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_P_Selct").release();
  Spectrum *sOtherNuNC_P_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_P_Selct").release();
  Spectrum *sCosmic_P_Selct = LoadFromFile<Spectrum>(fLoad,"sCosmic_P_Selct").release();
  Spectrum *sInTime_P_Selct = LoadFromFile<Spectrum>(fLoad,"sInTime_P_Selct").release();

  Spectrum *sSignal_TruP_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_TruP_Selct").release();
  Spectrum *sSignal_NotP_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_NotP_Selct").release();
  Spectrum *sBackgd_P_Selct = LoadFromFile<Spectrum>(fLoad,"sBackgd_P_Selct").release();

  Spectrum *sRecoVsTrue_ProtonP = LoadFromFile<Spectrum>(fLoad,"sRecoVsTrue_ProtonP").release();


  Spectrum *sAll_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sAll_PSimple_Selct").release();
  Spectrum *sSignal_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sSignal_PSimple_Selct").release();
  Spectrum *sOtherNuCC_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sOtherNuCC_PSimple_Selct").release();
  Spectrum *sOtherNuNC_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sOtherNuNC_PSimple_Selct").release();
  Spectrum *sCosmic_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sCosmic_PSimple_Selct").release();
  Spectrum *sInTime_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sInTime_PSimple_Selct").release();

  Spectrum *sSignal_TruPSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sSignal_TruPSimple_Selct").release();
  Spectrum *sSignal_NotPSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sSignal_NotPSimple_Selct").release();
  Spectrum *sBackgd_PSimple_Selct = LoadFromFile<Spectrum>(fLoad, "sBackgd_PSimple_Selct").release();

  Spectrum *sSignal_TruPSimple_Selct_Chi2P = LoadFromFile<Spectrum>(fLoad, "sSignal_TruPSimple_Selct_Chi2P").release();
  Spectrum *sSignal_TruPSimple_Selct_Chi2Mu = LoadFromFile<Spectrum>(fLoad, "sSignal_TruPSimple_Selct_Chi2Mu").release();
  Spectrum *sSignal_TruStopPSimple_Selct_Chi2P = LoadFromFile<Spectrum>(fLoad, "sSignal_TruStopPSimple_Selct_Chi2P").release();
  Spectrum *sSignal_TruStopPSimple_Selct_Chi2Mu = LoadFromFile<Spectrum>(fLoad, "sSignal_TruStopPSimple_Selct_Chi2Mu").release();



  Spectrum *sSignal_PMultPreFSI_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPreFSI_NoCut").release();
  Spectrum *sSignal_PMultPostFSI_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSI_NoCut").release();
  Spectrum *sSignal_PMultPostFSIMinE_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSIMinE_NoCut").release();
  Spectrum *sSignal_PMultReco_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultReco_NoCut").release();
  Spectrum *sSignal_PMultRecoStubs_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultRecoStubs_NoCut").release();

  Spectrum *sSignal_PMultPreFSI_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPreFSI_Selct").release();
  Spectrum *sSignal_PMultPostFSI_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSI_Selct").release();
  Spectrum *sSignal_PMultPostFSIMinE_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSIMinE_Selct").release();
  Spectrum *sSignal_PMultReco_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultReco_Selct").release();
  Spectrum *sSignal_PMultRecoStubs_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultRecoStubs_Selct").release();

  Spectrum *sRecoVsTrue_ProtonMult_NoCut = LoadFromFile<Spectrum>(fLoad,"sRecoVsTrue_ProtonMult_NoCut").release();
  Spectrum *sRecoVsTrue_ProtonMult_Selct = LoadFromFile<Spectrum>(fLoad,"sRecoVsTrue_ProtonMult_Selct").release();


  Spectrum *sSignal_PMultPostFSIE100_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSIE100_NoCut").release();
  Spectrum *sSignal_PMultPostFSIE100_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultPostFSIE100_Selct").release();
  Spectrum *sSignal_PMultRecoCheated_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultRecoCheated_NoCut").release();
  Spectrum *sSignal_PMultRecoCheated_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_PMultRecoCheated_Selct").release();

  Spectrum *sSignal_TrueProtonKE_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonKE_SpillMV").release();
  Spectrum *sSignal_TrueProtonKE_CheatedReco_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonKE_CheatedReco_SpillMV").release();
  Spectrum *sSignal_TrueProtonKE_RES_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonKE_RES_SpillMV").release();
  Spectrum *sSignal_TrueProtonKE_RES_CheatedReco_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonKE_RES_CheatedReco_SpillMV").release();

  Spectrum *sSignal_TrueProtonMom_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonMom_SpillMV").release();
  Spectrum *sSignal_TrueProtonMom_CheatedReco_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueProtonMom_CheatedReco_SpillMV").release();

  Spectrum *sSignal_TrueNuE_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueNuE_SpillMV").release();
  Spectrum *sSignal_TrueNuE_Selct_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueNuE_Selct_SpillMV").release();
  Spectrum *sSignalwProton_TrueMuMom_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_TrueMuMom_SpillMV").release();
  Spectrum *sSignalwProton_TrueMuMom_Selct_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_TrueMuMom_Selct_SpillMV").release();
  Spectrum *sSignalwProton_RecoMuMom_Selct_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_RecoMuMom_Selct_SpillMV").release();
  Spectrum *sAll_RecoMuMom_Selct_SpillMV = LoadFromFile<Spectrum>(fLoad,"sAll_RecoMuMom_Selct_SpillMV").release();
  Spectrum *sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic = LoadFromFile<Spectrum>(fLoad,"sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic").release();

  Spectrum *sSignal_TrueNuL_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueNuL_SpillMV").release();
  Spectrum *sSignal_TrueNuLOverE_SpillMV = LoadFromFile<Spectrum>(fLoad,"sSignal_TrueNuLOverE_SpillMV").release();

  Spectrum *sNumuCC_TrueNuE_SpillMV = LoadFromFile<Spectrum>(fLoad,"sNumuCC_TrueNuE_SpillMV").release();
  Spectrum *sNumuCC_TrueNuL_SpillMV = LoadFromFile<Spectrum>(fLoad,"sNumuCC_TrueNuL_SpillMV").release();
  Spectrum *sNumuCC_TrueNuLOverE_SpillMV = LoadFromFile<Spectrum>(fLoad,"sNumuCC_TrueNuLOverE_SpillMV").release();


  Spectrum *sNuMIMuons_CosThXZ = LoadFromFile<Spectrum>(fLoad,"sNuMIMuons_CosThXZ").release();
  Spectrum *sNuMIMuons_CosThNuMI = LoadFromFile<Spectrum>(fLoad,"sNuMIMuons_CosThNuMI").release();
  Spectrum *sNuMIMuons_CosThNuMINotIsoch = LoadFromFile<Spectrum>(fLoad,"sNuMIMuons_CosThNuMINotIsoch").release();
  Spectrum *sNuMIMuons_CosThNuMINotIsochNotPerp = LoadFromFile<Spectrum>(fLoad,"sNuMIMuons_CosThNuMINotIsochNotPerp").release();



  Spectrum *sAll_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sAll_Cont_Selct").release();
  Spectrum *sSignal_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sSignal_Cont_Selct").release();
  Spectrum *sOtherNuCC_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_Cont_Selct").release();
  Spectrum *sOtherNuNC_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_Cont_Selct").release();
  Spectrum *sCosmic_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sCosmic_Cont_Selct").release();
  Spectrum *sInTime_Cont_Selct = LoadFromFile<Spectrum>(fLoad,"sInTime_Cont_Selct").release();

  Spectrum *sAll_SelctG = LoadFromFile<Spectrum>(fLoad,"sAll_SelctG").release();
  Spectrum *sSignal_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignal_SelctG").release();
  Spectrum *sOtherNuCC_SelctG = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_SelctG").release();
  Spectrum *sOtherNuNC_SelctG = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_SelctG").release();
  Spectrum *sCosmic_SelctG = LoadFromFile<Spectrum>(fLoad,"sCosmic_SelctG").release();
  Spectrum *sInTime_SelctG = LoadFromFile<Spectrum>(fLoad,"sInTime_SelctG").release();

  Spectrum *sSignalQEL_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_SelctG").release();
  Spectrum *sSignalMEC_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_SelctG").release();
  Spectrum *sSignalRES_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalRES_SelctG").release();
  Spectrum *sSignalDIS_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_SelctG").release();
  Spectrum *sSignalCOH_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_SelctG").release();
  Spectrum *sSignalELS_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignalELS_SelctG").release();

  // GRAY's selection with proton vars
  Spectrum *sAll_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sAll_P_SelctG").release();
  Spectrum *sSignal_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignal_P_SelctG").release();
  Spectrum *sOtherNuCC_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_P_SelctG").release();
  Spectrum *sOtherNuNC_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_P_SelctG").release();
  Spectrum *sCosmic_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sCosmic_P_SelctG").release();
  Spectrum *sInTime_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sInTime_P_SelctG").release();

  Spectrum *sSignal_TruP_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignal_TruP_SelctG").release();
  Spectrum *sSignal_NotP_SelctG = LoadFromFile<Spectrum>(fLoad,"sSignal_NotP_SelctG").release();
  Spectrum *sBackgd_P_SelctG = LoadFromFile<Spectrum>(fLoad,"sBackgd_P_SelctG").release();


  Spectrum *sAll_SelTot = LoadFromFile<Spectrum>(fLoad,"sAll_SelTot").release();
  Spectrum *sSignal_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignal_SelTot").release();
  Spectrum *sOtherNuCC_SelTot = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_SelTot").release();
  Spectrum *sOtherNuNC_SelTot = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_SelTot").release();
  Spectrum *sCosmic_SelTot = LoadFromFile<Spectrum>(fLoad,"sCosmic_SelTot").release();
  Spectrum *sInTime_SelTot = LoadFromFile<Spectrum>(fLoad,"sInTime_SelTot").release();
  Spectrum *sSignalQEL_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_SelTot").release();
  Spectrum *sSignalMEC_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_SelTot").release();
  Spectrum *sSignalRES_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalRES_SelTot").release();
  Spectrum *sSignalDIS_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_SelTot").release();
  Spectrum *sSignalCOH_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_SelTot").release();
  Spectrum *sSignalELS_SelTot = LoadFromFile<Spectrum>(fLoad,"sSignalELS_SelTot").release();

  Spectrum *sSignalQEL_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_NoCut").release();
  Spectrum *sSignalMEC_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_NoCut").release();
  Spectrum *sSignalRES_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalRES_NoCut").release();
  Spectrum *sSignalDIS_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_NoCut").release();
  Spectrum *sSignalCOH_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_NoCut").release();
  //Spectrum *sSignalELS_NoCut = LoadFromFile<Spectrum>(fLoad,"sSignalELS_NoCut").release();

  Spectrum *sAll_SelCoh = LoadFromFile<Spectrum>(fLoad,"sAll_SelCoh").release();
  Spectrum *sSignal_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignal_SelCoh").release();
  Spectrum *sOtherNuCC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_SelCoh").release();
  Spectrum *sOtherNuNC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_SelCoh").release();
  Spectrum *sCosmic_SelCoh = LoadFromFile<Spectrum>(fLoad,"sCosmic_SelCoh").release();
  Spectrum *sInTime_SelCoh = LoadFromFile<Spectrum>(fLoad,"sInTime_SelCoh").release();
  Spectrum *sSignalQEL_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_SelCoh").release();
  Spectrum *sSignalMEC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_SelCoh").release();
  Spectrum *sSignalRES_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalRES_SelCoh").release();
  Spectrum *sSignalDIS_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_SelCoh").release();
  Spectrum *sSignalCOH_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_SelCoh").release();
  //Spectrum *sSignalELS_SelCoh = LoadFromFile<Spectrum>(fLoad,"sSignalELS_SelCoh").release();

  Spectrum *sT_All_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_All_SelCoh").release();
  Spectrum *sT_Signal_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_Signal_SelCoh").release();
  Spectrum *sT_OtherNuCC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_OtherNuCC_SelCoh").release();
  Spectrum *sT_OtherNuNC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_OtherNuNC_SelCoh").release();
  Spectrum *sT_Cosmic_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_Cosmic_SelCoh").release();
  Spectrum *sT_InTime_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_InTime_SelCoh").release();
  Spectrum *sT_SignalQEL_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_SignalQEL_SelCoh").release();
  Spectrum *sT_SignalMEC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_SignalMEC_SelCoh").release();
  Spectrum *sT_SignalRES_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_SignalRES_SelCoh").release();
  Spectrum *sT_SignalDIS_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_SignalDIS_SelCoh").release();
  Spectrum *sT_SignalCOH_SelCoh = LoadFromFile<Spectrum>(fLoad,"sT_SignalCOH_SelCoh").release();

  Spectrum *sTzoom_All_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_All_SelCoh").release();
  Spectrum *sTzoom_Signal_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Signal_SelCoh").release();
  Spectrum *sTzoom_OtherNuCC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuCC_SelCoh").release();
  Spectrum *sTzoom_OtherNuNC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuNC_SelCoh").release();
  Spectrum *sTzoom_Cosmic_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Cosmic_SelCoh").release();
  Spectrum *sTzoom_InTime_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_InTime_SelCoh").release();
  Spectrum *sTzoom_SignalQEL_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalQEL_SelCoh").release();
  Spectrum *sTzoom_SignalMEC_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalMEC_SelCoh").release();
  Spectrum *sTzoom_SignalRES_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalRES_SelCoh").release();
  Spectrum *sTzoom_SignalDIS_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalDIS_SelCoh").release();
  Spectrum *sTzoom_SignalCOH_SelCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalCOH_SelCoh").release();

  Spectrum *sTzoom_All_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_All_ChtCoh").release();
  Spectrum *sTzoom_Signal_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Signal_ChtCoh").release();
  Spectrum *sTzoom_OtherNuCC_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuCC_ChtCoh").release();
  Spectrum *sTzoom_OtherNuNC_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuNC_ChtCoh").release();
  Spectrum *sTzoom_Cosmic_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Cosmic_ChtCoh").release();
  Spectrum *sTzoom_InTime_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_InTime_ChtCoh").release();
  Spectrum *sTzoom_SignalQEL_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalQEL_ChtCoh").release();
  Spectrum *sTzoom_SignalMEC_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalMEC_ChtCoh").release();
  Spectrum *sTzoom_SignalRES_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalRES_ChtCoh").release();
  Spectrum *sTzoom_SignalDIS_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalDIS_ChtCoh").release();
  Spectrum *sTzoom_SignalCOH_ChtCoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalCOH_ChtCoh").release();
  //
  Spectrum *sTzoom_All_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_All_CLECoh").release();
  Spectrum *sTzoom_Signal_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Signal_CLECoh").release();
  Spectrum *sTzoom_OtherNuCC_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuCC_CLECoh").release();
  Spectrum *sTzoom_OtherNuNC_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_OtherNuNC_CLECoh").release();
  Spectrum *sTzoom_Cosmic_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_Cosmic_CLECoh").release();
  Spectrum *sTzoom_InTime_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_InTime_CLECoh").release();
  Spectrum *sTzoom_SignalQEL_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalQEL_CLECoh").release();
  Spectrum *sTzoom_SignalMEC_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalMEC_CLECoh").release();
  Spectrum *sTzoom_SignalRES_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalRES_CLECoh").release();
  Spectrum *sTzoom_SignalDIS_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalDIS_CLECoh").release();
  Spectrum *sTzoom_SignalCOH_CLECoh = LoadFromFile<Spectrum>(fLoad,"sTzoom_SignalCOH_CLECoh").release();

  Spectrum *sAll_SelRes = LoadFromFile<Spectrum>(fLoad,"sAll_SelRes").release();
  Spectrum *sSignal_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignal_SelRes").release();
  Spectrum *sOtherNuCC_SelRes = LoadFromFile<Spectrum>(fLoad,"sOtherNuCC_SelRes").release();
  Spectrum *sOtherNuNC_SelRes = LoadFromFile<Spectrum>(fLoad,"sOtherNuNC_SelRes").release();
  Spectrum *sCosmic_SelRes = LoadFromFile<Spectrum>(fLoad,"sCosmic_SelRes").release();
  Spectrum *sInTime_SelRes = LoadFromFile<Spectrum>(fLoad,"sInTime_SelRes").release();
  Spectrum *sSignalQEL_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignalQEL_SelRes").release();
  Spectrum *sSignalMEC_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignalMEC_SelRes").release();
  Spectrum *sSignalRES_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignalRES_SelRes").release();
  Spectrum *sSignalDIS_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignalDIS_SelRes").release();
  Spectrum *sSignalCOH_SelRes = LoadFromFile<Spectrum>(fLoad,"sSignalCOH_SelRes").release();

  // Spectra for 1MuNP selection
  Spectrum *sSignalwProton_recoMuMom_BaseSliceSel = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_BaseSliceSel").release();
  Spectrum *sSignalwProton_recoMuMom_TrueSignal = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_TrueSignal").release();
  Spectrum *sSignalwProton_recoMuMom_NumuMuWrong = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NumuMuWrong").release();
  Spectrum *sSignalwProton_recoMuMom_NumuPWrong = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NumuPWrong").release();
  Spectrum *sSignalwProton_recoMuMom_NumuBothWrong = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NumuBothWrong").release();
  Spectrum *sSignalwProton_recoMuMom_NumuNonFidVol = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NumuNonFidVol").release();
  Spectrum *sSignalwProton_recoMuMom_NumuOther = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NumuOther").release();
  Spectrum *sSignalwProton_recoMuMom_NuOther = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_NuOther").release();
  Spectrum *sSignalwProton_recoMuMom_Cosmic = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_Cosmic").release();
  Spectrum *sSignalwProton_recoMuMom_Cosmic_InTime = LoadFromFile<Spectrum>(fLoad,"sSignalwProton_recoMuMom_Cosmic_InTime").release();


  // Scale in-time cosmic background (see https://github.com/SBNSoftware/sbnana/blob/develop/sbnana/SBNAna/anademo/demo_exposure.C)
  const double nomIntensity = 5.0e13;
  const double nomLivetime = pot / nomIntensity;
  const double beamLivetime = sSignal_NoCut->FakeData(pot).Livetime();
  const double cosmicLivetime = nomLivetime - beamLivetime;
  /////////////

  TH1* hAll_NoCut      = sAll_NoCut->ToTH1( pot, kBlack );
  hAll_NoCut->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hAll_Selct      = sAll_Selct->ToTH1( pot, kBlack );
  hAll_Selct->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hAll_P_Selct      = sAll_P_Selct->ToTH1( pot, kBlack );
  hAll_P_Selct->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hAll_PSimple_Selct= sAll_PSimple_Selct->ToTH1( pot, kBlack );
  hAll_PSimple_Selct->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hAll_Cont_Selct = sAll_Cont_Selct->ToTH1( pot, kBlack );
  hAll_Cont_Selct->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hAll_SelctG     = sAll_SelctG->ToTH1( pot, kBlack );
  hAll_SelctG->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );


  TH1* hAll_P_SelctG      = sAll_P_SelctG->ToTH1( pot, kBlack );
  hAll_P_SelctG->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );


  TH1* hAll_SelTot     = sAll_SelTot->ToTH1( pot, kBlack );
  hAll_SelTot->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hAll_SelCoh     = sAll_SelCoh->ToTH1( pot, kBlack );
  hAll_SelCoh->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hT_All_SelCoh     = sT_All_SelCoh->ToTH1( pot, kBlack );
  hT_All_SelCoh->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hTzoom_All_SelCoh     = sTzoom_All_SelCoh->ToTH1( pot, kBlack );
  hTzoom_All_SelCoh->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hTzoom_All_ChtCoh     = sTzoom_All_ChtCoh->ToTH1( pot, kBlack );
  hTzoom_All_ChtCoh->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );
  TH1* hTzoom_All_CLECoh     = sTzoom_All_CLECoh->ToTH1( pot, kBlack );
  hTzoom_All_CLECoh->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hAll_SelRes     = sAll_SelRes->ToTH1( pot, kBlack );
  hAll_SelRes->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );

  TH1* hSignalwProton_recoMuMom_BaseSliceSel = sSignalwProton_recoMuMom_BaseSliceSel->ToTH1( pot, kBlack );
  hSignalwProton_recoMuMom_BaseSliceSel->GetYaxis()->SetTitle( TString::Format("Slices / %.3fe20 POT",pot/1.0e20) );


  TH1* hSignal_NoCut    = sSignal_NoCut->ToTH1( pot, colorwheel[0] );    hSignal_NoCut->SetFillColor(colorwheel[0]);
  TH1* hSignal_Selct    = sSignal_Selct->ToTH1( pot, colorwheel[0] );    hSignal_Selct->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_NoCut = sOtherNuCC_NoCut->ToTH1( pot, colorwheel[1] ); hOtherNuCC_NoCut->SetFillColor(colorwheel[1]);
  TH1* hOtherNuCC_Selct = sOtherNuCC_Selct->ToTH1( pot, colorwheel[1] ); hOtherNuCC_Selct->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_NoCut = sOtherNuNC_NoCut->ToTH1( pot, colorwheel[2] ); hOtherNuNC_NoCut->SetFillColor(colorwheel[2]);
  TH1* hOtherNuNC_Selct = sOtherNuNC_Selct->ToTH1( pot, colorwheel[2] ); hOtherNuNC_Selct->SetFillColor(colorwheel[2]);
  TH1* hCosmic_NoCut    = sCosmic_NoCut->ToTH1( pot, colorwheel[3] );    hCosmic_NoCut->SetFillColor(colorwheel[3]);
  TH1* hCosmic_Selct    = sCosmic_Selct->ToTH1( pot, colorwheel[3] );    hCosmic_Selct->SetFillColor(colorwheel[3]);

  TH1* hSignal_P_Selct    = sSignal_P_Selct->ToTH1( pot, colorwheel[0] );    hSignal_P_Selct->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_P_Selct = sOtherNuCC_P_Selct->ToTH1( pot, colorwheel[1] ); hOtherNuCC_P_Selct->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_P_Selct = sOtherNuNC_P_Selct->ToTH1( pot, colorwheel[2] ); hOtherNuNC_P_Selct->SetFillColor(colorwheel[2]);
  TH1* hCosmic_P_Selct    = sCosmic_P_Selct->ToTH1( pot, colorwheel[3] );    hCosmic_P_Selct->SetFillColor(colorwheel[3]);
  TH1* hInTime_P_Selct    = sInTime_P_Selct->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_P_Selct->SetFillColor(colorwheel[4]);

  TH1* hSignal_TruP_Selct = sSignal_TruP_Selct->ToTH1( pot, colorwheel[0] ); hSignal_TruP_Selct->SetFillColor(colorwheel[0]);
  TH1* hSignal_NotP_Selct = sSignal_NotP_Selct->ToTH1( pot, colorwheel[0] ); hSignal_NotP_Selct->SetFillColor(colorwheel[1]);
  TH1* hBackgd_P_Selct    = sBackgd_P_Selct->ToTH1( pot, colorwheel[0] );    hBackgd_P_Selct->SetFillColor(colorwheel[2]);

  TH2* hRecoVsTrue_ProtonP = sRecoVsTrue_ProtonP->ToTH2( pot );


  TH1* hSignal_PSimple_Selct    = sSignal_PSimple_Selct->ToTH1( pot, colorwheel[0] );    hSignal_PSimple_Selct->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_PSimple_Selct = sOtherNuCC_PSimple_Selct->ToTH1( pot, colorwheel[1] ); hOtherNuCC_PSimple_Selct->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_PSimple_Selct = sOtherNuNC_PSimple_Selct->ToTH1( pot, colorwheel[2] ); hOtherNuNC_PSimple_Selct->SetFillColor(colorwheel[2]);
  TH1* hCosmic_PSimple_Selct    = sCosmic_PSimple_Selct->ToTH1( pot, colorwheel[3] );    hCosmic_PSimple_Selct->SetFillColor(colorwheel[3]);
  TH1* hInTime_PSimple_Selct    = sInTime_PSimple_Selct->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_PSimple_Selct->SetFillColor(colorwheel[4]);

  TH1* hSignal_TruPSimple_Selct = sSignal_TruPSimple_Selct->ToTH1( pot, colorwheel[0] ); hSignal_TruPSimple_Selct->SetFillColor(colorwheel[0]);
  TH1* hSignal_NotPSimple_Selct = sSignal_NotPSimple_Selct->ToTH1( pot, colorwheel[1] ); hSignal_NotPSimple_Selct->SetFillColor(colorwheel[1]);
  TH1* hBackgd_PSimple_Selct    = sBackgd_PSimple_Selct->ToTH1( pot, colorwheel[2] );    hBackgd_PSimple_Selct->SetFillColor(colorwheel[2]);

  TH1* hSignal_TruPSimple_Selct_Chi2P = sSignal_TruPSimple_Selct_Chi2P->ToTH1( pot, kBlack );
  TH1* hSignal_TruPSimple_Selct_Chi2Mu = sSignal_TruPSimple_Selct_Chi2Mu->ToTH1( pot, kBlack );
  TH1* hSignal_TruStopPSimple_Selct_Chi2P = sSignal_TruStopPSimple_Selct_Chi2P->ToTH1( pot, kBlack );
  TH1* hSignal_TruStopPSimple_Selct_Chi2Mu = sSignal_TruStopPSimple_Selct_Chi2Mu->ToTH1( pot, kBlack );



  // Proton Multiplicity
  // PLOT 1 - Truth info
  // No Cut
  TH1* hSignal_PMultPreFSI_NoCut = sSignal_PMultPreFSI_NoCut->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultPostFSI_NoCut = sSignal_PMultPostFSI_NoCut->ToTH1( pot, colorwheel_mode[1] );
  TH1* hSignal_PMultPostFSIMinE_NoCut = sSignal_PMultPostFSIMinE_NoCut->ToTH1( pot, colorwheel_mode[2] );
  TH1* hSignal_PMultPostFSIE100_NoCut = sSignal_PMultPostFSIE100_NoCut->ToTH1( pot, colorwheel_mode[3] );
  // Selct
  TH1* hSignal_PMultPreFSI_Selct = sSignal_PMultPreFSI_Selct->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultPostFSI_Selct = sSignal_PMultPostFSI_Selct->ToTH1( pot, colorwheel_mode[1] );
  TH1* hSignal_PMultPostFSIMinE_Selct = sSignal_PMultPostFSIMinE_Selct->ToTH1( pot, colorwheel_mode[2] );
  TH1* hSignal_PMultPostFSIE100_Selct = sSignal_PMultPostFSIE100_Selct->ToTH1( pot, colorwheel_mode[3] );

  // PLOT 2 - Reco w/o Stubs, compare to 100MeV cut... ?
  // No Cut
  TH1* hSignal_PMultPostFSIE100_RecoComp_NoCut = sSignal_PMultPostFSIE100_NoCut->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultReco_NoCut = sSignal_PMultReco_NoCut->ToTH1( pot, colorwheel_mode[1] );
  TH1* hSignal_PMultRecoCheated_NoCut = sSignal_PMultRecoCheated_NoCut->ToTH1( pot, colorwheel_mode[2] );
  // Selct
  TH1* hSignal_PMultPostFSIE100_RecoComp_Selct = sSignal_PMultPostFSIE100_Selct->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultReco_Selct = sSignal_PMultReco_Selct->ToTH1( pot, colorwheel_mode[1] );
  TH1* hSignal_PMultRecoCheated_Selct = sSignal_PMultRecoCheated_Selct->ToTH1( pot, colorwheel_mode[2] );

  // PLOT 3 - Reco w/ stubs, compare to 15 MeV version... ?
  // No Cut
  TH1* hSignal_PMultPostFSIMinE_RecoComp_NoCut = sSignal_PMultPostFSIMinE_NoCut->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultRecoStubs_NoCut = sSignal_PMultRecoStubs_NoCut->ToTH1( pot, colorwheel_mode[1] );
  // Selct
  TH1* hSignal_PMultPostFSIMinE_RecoComp_Selct = sSignal_PMultPostFSIMinE_Selct->ToTH1( pot, colorwheel_mode[0] );
  TH1* hSignal_PMultRecoStubs_Selct = sSignal_PMultRecoStubs_Selct->ToTH1( pot, colorwheel_mode[1] );


  TH2* hRecoVsTrue_ProtonMult_NoCut = sRecoVsTrue_ProtonMult_NoCut->ToTH2( pot );
  TH2* hRecoVsTrue_ProtonMult_Selct = sRecoVsTrue_ProtonMult_Selct->ToTH2( pot );

  // True multiplicity and efficiency, as energy calculation...
  TH1* hSignal_TrueProtonKE_SpillMV = sSignal_TrueProtonKE_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignal_TrueProtonKE_CheatedReco_SpillMV = sSignal_TrueProtonKE_CheatedReco_SpillMV->ToTH1( pot, kRed );
  Ratio rsSignal_TrueProtonKE_SpillMV( *sSignal_TrueProtonKE_CheatedReco_SpillMV, *sSignal_TrueProtonKE_SpillMV, true );
  TH1* rSignal_TrueProtonKE_SpillMV = rsSignal_TrueProtonKE_SpillMV.ToTH1( kBlack );
 
  TH1* hSignal_TrueProtonKE_RES_SpillMV = sSignal_TrueProtonKE_RES_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignal_TrueProtonKE_RES_CheatedReco_SpillMV = sSignal_TrueProtonKE_RES_CheatedReco_SpillMV->ToTH1( pot, kRed );
  Ratio rsSignal_TrueProtonKE_RES_SpillMV( *sSignal_TrueProtonKE_RES_CheatedReco_SpillMV, *sSignal_TrueProtonKE_RES_SpillMV, true );
  TH1* rSignal_TrueProtonKE_RES_SpillMV = rsSignal_TrueProtonKE_RES_SpillMV.ToTH1( kBlack );

  TH1* hSignal_TrueProtonMom_SpillMV = sSignal_TrueProtonMom_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignal_TrueProtonMom_CheatedReco_SpillMV = sSignal_TrueProtonMom_CheatedReco_SpillMV->ToTH1( pot, kRed );
  Ratio rsSignal_TrueProtonMom_SpillMV( *sSignal_TrueProtonMom_CheatedReco_SpillMV, *sSignal_TrueProtonMom_SpillMV, true );
  TH1* rSignal_TrueProtonMom_SpillMV = rsSignal_TrueProtonMom_SpillMV.ToTH1( kBlack );

  TH1* hSignal_TrueNuE_SpillMV = sSignal_TrueNuE_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignal_TrueNuE_Selct_SpillMV = sSignal_TrueNuE_Selct_SpillMV->ToTH1( pot, kRed );
  Ratio rsSignal_TrueNuE_SpillMV( *sSignal_TrueNuE_Selct_SpillMV, *sSignal_TrueNuE_SpillMV, true );
  TH1* rSignal_TrueNuE_SpillMV = rsSignal_TrueNuE_SpillMV.ToTH1( kBlack );

  TH1* hSignalwProton_TrueMuMom_SpillMV = sSignalwProton_TrueMuMom_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignalwProton_TrueMuMom_Selct_SpillMV = sSignalwProton_TrueMuMom_Selct_SpillMV->ToTH1( pot, kRed );
  Ratio rsSignalwProton_TrueMuMom_SpillMV( *sSignalwProton_TrueMuMom_Selct_SpillMV, *sSignalwProton_TrueMuMom_SpillMV, true );
  TH1* rSignalwProton_TrueMuMom_SpillMV = rsSignalwProton_TrueMuMom_SpillMV.ToTH1( kBlack );

  TH1* hSignalwProton_RecoMuMom_Selct_SpillMV = sSignalwProton_RecoMuMom_Selct_SpillMV->ToTH1( pot, kBlack );
  TH1* rSignalwProton_RecoMuMom_Selct_SpillMV = sSignalwProton_RecoMuMom_Selct_SpillMV->ToTH1( pot, kBlack );
  TH1* hAll_RecoMuMom_Selct_SpillMV = sAll_RecoMuMom_Selct_SpillMV->ToTH1( pot, kRed );
  hAll_RecoMuMom_Selct_SpillMV->Add( sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic->ToTH1(cosmicLivetime, kRed, kSolid, kLivetime ) );
  // Here using Ratio doesn't quite work because the exposure definitions are different...
  //   Ratio rsSignalwProton_RecoMuMom_Selct_SpillMV( *sSignalwProton_RecoMuMom_Selct_SpillMV, *sAll_RecoMuMom_Selct_SpillMV + *sAll_RecoMuMom_Selct_SpillMV_InTimeCosmic, true );
  //   TH1* rSignalwProton_RecoMuMom_Selct_SpillMV = rsSignalwProton_RecoMuMom_Selct_SpillMV.ToTH1( kRed );
  // Instead, we will divide histograms
  rSignalwProton_RecoMuMom_Selct_SpillMV->Divide(rSignalwProton_RecoMuMom_Selct_SpillMV,hAll_RecoMuMom_Selct_SpillMV,1,1,"B");

  TH1* hSignal_TrueNuL_SpillMV = sSignal_TrueNuL_SpillMV->ToTH1( pot, kBlack );
  TH1* hSignal_TrueNuLOverE_SpillMV = sSignal_TrueNuLOverE_SpillMV->ToTH1( pot, kBlack );

  TH1* hNumuCC_TrueNuE_SpillMV = sNumuCC_TrueNuE_SpillMV->ToTH1( pot, kBlack );
  TH1* hNumuCC_TrueNuL_SpillMV = sNumuCC_TrueNuL_SpillMV->ToTH1( pot, kBlack );
  TH1* hNumuCC_TrueNuLOverE_SpillMV = sNumuCC_TrueNuLOverE_SpillMV->ToTH1( pot, kBlack );


  TH1* hNuMIMuons_CosThXZ = sNuMIMuons_CosThXZ->ToTH1( pot, kBlack );
  TH1* hNuMIMuons_CosThNuMI = sNuMIMuons_CosThNuMI->ToTH1( pot, kBlack );
  TH1* hNuMIMuons_CosThNuMINotIsoch = sNuMIMuons_CosThNuMINotIsoch->ToTH1( pot, kRed );
  TH1* hNuMIMuons_CosThNuMINotIsochNotPerp = sNuMIMuons_CosThNuMINotIsoch->ToTH1( pot, kGreen+2 );


  TH1* hSignalQEL_Selct = sSignalQEL_Selct->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_Selct->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_Selct = sSignalMEC_Selct->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_Selct->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_Selct = sSignalRES_Selct->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_Selct->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_Selct = sSignalDIS_Selct->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_Selct->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_Selct = sSignalCOH_Selct->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_Selct->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_Selct = sOtherNuCC_Selct->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_Selct->SetFillColor(colorwheel_mode[5]);

  TH1* hInTime_NoCut    = sInTime_NoCut->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_NoCut->SetFillColor(colorwheel[4]);
  TH1* hInTime_Selct    = sInTime_Selct->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_Selct->SetFillColor(colorwheel[4]);

  TH1* hSignal_Cont_NoCut    = sSignal_Cont_NoCut->ToTH1( pot, colorwheel[0] );    hSignal_Cont_NoCut->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_Cont_NoCut = sOtherNuCC_Cont_NoCut->ToTH1( pot, colorwheel[1] ); hOtherNuCC_Cont_NoCut->SetFillColor(colorwheel[1]);

  TH1* hSignal_Cont_Selct    = sSignal_Cont_Selct->ToTH1( pot, colorwheel[0] );    hSignal_Cont_Selct->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_Cont_Selct = sOtherNuCC_Cont_Selct->ToTH1( pot, colorwheel[1] ); hOtherNuCC_Cont_Selct->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_Cont_Selct = sOtherNuNC_Cont_Selct->ToTH1( pot, colorwheel[2] ); hOtherNuNC_Cont_Selct->SetFillColor(colorwheel[2]);
  TH1* hCosmic_Cont_Selct    = sCosmic_Cont_Selct->ToTH1( pot, colorwheel[3] );    hCosmic_Cont_Selct->SetFillColor(colorwheel[3]);
  TH1* hInTime_Cont_Selct    = sInTime_Cont_Selct->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_Cont_Selct->SetFillColor(colorwheel[4]);

  TH1* hSignal_SelctG    = sSignal_SelctG->ToTH1( pot, colorwheel[0] );    hSignal_SelctG->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_SelctG = sOtherNuCC_SelctG->ToTH1( pot, colorwheel[1] ); hOtherNuCC_SelctG->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_SelctG = sOtherNuNC_SelctG->ToTH1( pot, colorwheel[2] ); hOtherNuNC_SelctG->SetFillColor(colorwheel[2]);
  TH1* hCosmic_SelctG    = sCosmic_SelctG->ToTH1( pot, colorwheel[3] );    hCosmic_SelctG->SetFillColor(colorwheel[3]);
  TH1* hInTime_SelctG    = sInTime_SelctG->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_SelctG->SetFillColor(colorwheel[4]);

  TH1* hSignalQEL_SelctG = sSignalQEL_SelctG->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_SelctG->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_SelctG = sSignalMEC_SelctG->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_SelctG->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_SelctG = sSignalRES_SelctG->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_SelctG->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_SelctG = sSignalDIS_SelctG->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_SelctG->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_SelctG = sSignalCOH_SelctG->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_SelctG->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_SelctG = sOtherNuCC_SelctG->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_SelctG->SetFillColor(colorwheel_mode[5]);


  TH1* hSignal_P_SelctG    = sSignal_P_SelctG->ToTH1( pot, colorwheel[0] );    hSignal_P_SelctG->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_P_SelctG = sOtherNuCC_P_SelctG->ToTH1( pot, colorwheel[1] ); hOtherNuCC_P_SelctG->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_P_SelctG = sOtherNuNC_P_SelctG->ToTH1( pot, colorwheel[2] ); hOtherNuNC_P_SelctG->SetFillColor(colorwheel[2]);
  TH1* hCosmic_P_SelctG    = sCosmic_P_SelctG->ToTH1( pot, colorwheel[3] );    hCosmic_P_SelctG->SetFillColor(colorwheel[3]);
  TH1* hInTime_P_SelctG    = sInTime_P_SelctG->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_P_SelctG->SetFillColor(colorwheel[4]);

  TH1* hSignal_TruP_SelctG = sSignal_TruP_SelctG->ToTH1( pot, colorwheel[0] ); hSignal_TruP_SelctG->SetFillColor(colorwheel[0]);
  TH1* hSignal_NotP_SelctG = sSignal_NotP_SelctG->ToTH1( pot, colorwheel[1] ); hSignal_NotP_SelctG->SetFillColor(colorwheel[1]);
  TH1* hBackgd_P_SelctG    = sBackgd_P_SelctG->ToTH1( pot, colorwheel[2] );    hBackgd_P_SelctG->SetFillColor(colorwheel[2]);


  TH1* hSignal_SelTot    = sSignal_SelTot->ToTH1( pot, colorwheel[0] );         hSignal_SelTot->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_SelTot = sOtherNuCC_SelTot->ToTH1( pot, colorwheel[1] );      hOtherNuCC_SelTot->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_SelTot = sOtherNuNC_SelTot->ToTH1( pot, colorwheel[2] );      hOtherNuNC_SelTot->SetFillColor(colorwheel[2]);
  TH1* hCosmic_SelTot    = sCosmic_SelTot->ToTH1( pot, colorwheel[3] );         hCosmic_SelTot->SetFillColor(colorwheel[3]);
  TH1* hInTime_SelTot    = sInTime_SelTot->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_SelTot->SetFillColor(colorwheel[4]);
  TH1* hSignalQEL_SelTot = sSignalQEL_SelTot->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_SelTot->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_SelTot = sSignalMEC_SelTot->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_SelTot->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_SelTot = sSignalRES_SelTot->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_SelTot->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_SelTot = sSignalDIS_SelTot->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_SelTot->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_SelTot = sSignalCOH_SelTot->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_SelTot->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_SelTot = sOtherNuCC_SelTot->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_SelTot->SetFillColor(colorwheel_mode[5]);

  TH1* hSignalQEL_NoCut = sSignalQEL_NoCut->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_NoCut->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_NoCut = sSignalMEC_NoCut->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_NoCut->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_NoCut = sSignalRES_NoCut->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_NoCut->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_NoCut = sSignalDIS_NoCut->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_NoCut->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_NoCut = sSignalCOH_NoCut->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_NoCut->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_NoCut = sOtherNuCC_NoCut->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_NoCut->SetFillColor(colorwheel_mode[5]);

  TH1* hSignal_SelCoh    = sSignal_SelCoh->ToTH1( pot, colorwheel[0] );         hSignal_SelCoh->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_SelCoh = sOtherNuCC_SelCoh->ToTH1( pot, colorwheel[1] );      hOtherNuCC_SelCoh->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_SelCoh = sOtherNuNC_SelCoh->ToTH1( pot, colorwheel[2] );      hOtherNuNC_SelCoh->SetFillColor(colorwheel[2]);
  TH1* hCosmic_SelCoh    = sCosmic_SelCoh->ToTH1( pot, colorwheel[3] );         hCosmic_SelCoh->SetFillColor(colorwheel[3]);
  TH1* hInTime_SelCoh    = sInTime_SelCoh->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_SelCoh->SetFillColor(colorwheel[4]);
  TH1* hSignalQEL_SelCoh = sSignalQEL_SelCoh->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_SelCoh->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_SelCoh = sSignalMEC_SelCoh->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_SelCoh->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_SelCoh = sSignalRES_SelCoh->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_SelCoh->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_SelCoh = sSignalDIS_SelCoh->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_SelCoh->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_SelCoh = sSignalCOH_SelCoh->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_SelCoh->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_SelCoh = sOtherNuCC_SelCoh->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_SelCoh->SetFillColor(colorwheel_mode[5]);

  TH1* hT_Signal_SelCoh    = sT_Signal_SelCoh->ToTH1( pot, colorwheel[0] );         hT_Signal_SelCoh->SetFillColor(colorwheel[0]);
  TH1* hT_OtherNuCC_SelCoh = sT_OtherNuCC_SelCoh->ToTH1( pot, colorwheel[1] );      hT_OtherNuCC_SelCoh->SetFillColor(colorwheel[1]);
  TH1* hT_OtherNuNC_SelCoh = sT_OtherNuNC_SelCoh->ToTH1( pot, colorwheel[2] );      hT_OtherNuNC_SelCoh->SetFillColor(colorwheel[2]);
  TH1* hT_Cosmic_SelCoh    = sT_Cosmic_SelCoh->ToTH1( pot, colorwheel[3] );         hT_Cosmic_SelCoh->SetFillColor(colorwheel[3]);
  TH1* hT_InTime_SelCoh    = sT_InTime_SelCoh->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hT_InTime_SelCoh->SetFillColor(colorwheel[4]);
  TH1* hT_SignalQEL_SelCoh = sT_SignalQEL_SelCoh->ToTH1( pot, colorwheel_mode[0] ); hT_SignalQEL_SelCoh->SetFillColor(colorwheel_mode[0]);
  TH1* hT_SignalMEC_SelCoh = sT_SignalMEC_SelCoh->ToTH1( pot, colorwheel_mode[1] ); hT_SignalMEC_SelCoh->SetFillColor(colorwheel_mode[1]);
  TH1* hT_SignalRES_SelCoh = sT_SignalRES_SelCoh->ToTH1( pot, colorwheel_mode[2] ); hT_SignalRES_SelCoh->SetFillColor(colorwheel_mode[2]);
  TH1* hT_SignalDIS_SelCoh = sT_SignalDIS_SelCoh->ToTH1( pot, colorwheel_mode[3] ); hT_SignalDIS_SelCoh->SetFillColor(colorwheel_mode[3]);
  TH1* hT_SignalCOH_SelCoh = sT_SignalCOH_SelCoh->ToTH1( pot, colorwheel_mode[4] ); hT_SignalCOH_SelCoh->SetFillColor(colorwheel_mode[4]);
  TH1* hT_SignalELS_SelCoh = sT_OtherNuCC_SelCoh->ToTH1( pot, colorwheel_mode[5] ); hT_SignalELS_SelCoh->SetFillColor(colorwheel_mode[5]);

  TH1* hTzoom_Signal_SelCoh    = sTzoom_Signal_SelCoh->ToTH1( pot, colorwheel[0] );         hTzoom_Signal_SelCoh->SetFillColor(colorwheel[0]);
  TH1* hTzoom_OtherNuCC_SelCoh = sTzoom_OtherNuCC_SelCoh->ToTH1( pot, colorwheel[1] );      hTzoom_OtherNuCC_SelCoh->SetFillColor(colorwheel[1]);
  TH1* hTzoom_OtherNuNC_SelCoh = sTzoom_OtherNuNC_SelCoh->ToTH1( pot, colorwheel[2] );      hTzoom_OtherNuNC_SelCoh->SetFillColor(colorwheel[2]);
  TH1* hTzoom_Cosmic_SelCoh    = sTzoom_Cosmic_SelCoh->ToTH1( pot, colorwheel[3] );         hTzoom_Cosmic_SelCoh->SetFillColor(colorwheel[3]);
  TH1* hTzoom_InTime_SelCoh    = sTzoom_InTime_SelCoh->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hTzoom_InTime_SelCoh->SetFillColor(colorwheel[4]);
  TH1* hTzoom_SignalQEL_SelCoh = sTzoom_SignalQEL_SelCoh->ToTH1( pot, colorwheel_mode[0] ); hTzoom_SignalQEL_SelCoh->SetFillColor(colorwheel_mode[0]);
  TH1* hTzoom_SignalMEC_SelCoh = sTzoom_SignalMEC_SelCoh->ToTH1( pot, colorwheel_mode[1] ); hTzoom_SignalMEC_SelCoh->SetFillColor(colorwheel_mode[1]);
  TH1* hTzoom_SignalRES_SelCoh = sTzoom_SignalRES_SelCoh->ToTH1( pot, colorwheel_mode[2] ); hTzoom_SignalRES_SelCoh->SetFillColor(colorwheel_mode[2]);
  TH1* hTzoom_SignalDIS_SelCoh = sTzoom_SignalDIS_SelCoh->ToTH1( pot, colorwheel_mode[3] ); hTzoom_SignalDIS_SelCoh->SetFillColor(colorwheel_mode[3]);
  TH1* hTzoom_SignalCOH_SelCoh = sTzoom_SignalCOH_SelCoh->ToTH1( pot, colorwheel_mode[4] ); hTzoom_SignalCOH_SelCoh->SetFillColor(colorwheel_mode[4]);
  TH1* hTzoom_SignalELS_SelCoh = sTzoom_OtherNuCC_SelCoh->ToTH1( pot, colorwheel_mode[5] ); hTzoom_SignalELS_SelCoh->SetFillColor(colorwheel_mode[5]);

  TH1* hTzoom_Signal_ChtCoh    = sTzoom_Signal_ChtCoh->ToTH1( pot, colorwheel[0] );         hTzoom_Signal_ChtCoh->SetFillColor(colorwheel[0]);
  TH1* hTzoom_OtherNuCC_ChtCoh = sTzoom_OtherNuCC_ChtCoh->ToTH1( pot, colorwheel[1] );      hTzoom_OtherNuCC_ChtCoh->SetFillColor(colorwheel[1]);
  TH1* hTzoom_OtherNuNC_ChtCoh = sTzoom_OtherNuNC_ChtCoh->ToTH1( pot, colorwheel[2] );      hTzoom_OtherNuNC_ChtCoh->SetFillColor(colorwheel[2]);
  TH1* hTzoom_Cosmic_ChtCoh    = sTzoom_Cosmic_ChtCoh->ToTH1( pot, colorwheel[3] );         hTzoom_Cosmic_ChtCoh->SetFillColor(colorwheel[3]);
  TH1* hTzoom_InTime_ChtCoh    = sTzoom_InTime_ChtCoh->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hTzoom_InTime_ChtCoh->SetFillColor(colorwheel[4]);
  TH1* hTzoom_SignalQEL_ChtCoh = sTzoom_SignalQEL_ChtCoh->ToTH1( pot, colorwheel_mode[0] ); hTzoom_SignalQEL_ChtCoh->SetFillColor(colorwheel_mode[0]);
  TH1* hTzoom_SignalMEC_ChtCoh = sTzoom_SignalMEC_ChtCoh->ToTH1( pot, colorwheel_mode[1] ); hTzoom_SignalMEC_ChtCoh->SetFillColor(colorwheel_mode[1]);
  TH1* hTzoom_SignalRES_ChtCoh = sTzoom_SignalRES_ChtCoh->ToTH1( pot, colorwheel_mode[2] ); hTzoom_SignalRES_ChtCoh->SetFillColor(colorwheel_mode[2]);
  TH1* hTzoom_SignalDIS_ChtCoh = sTzoom_SignalDIS_ChtCoh->ToTH1( pot, colorwheel_mode[3] ); hTzoom_SignalDIS_ChtCoh->SetFillColor(colorwheel_mode[3]);
  TH1* hTzoom_SignalCOH_ChtCoh = sTzoom_SignalCOH_ChtCoh->ToTH1( pot, colorwheel_mode[4] ); hTzoom_SignalCOH_ChtCoh->SetFillColor(colorwheel_mode[4]);
  TH1* hTzoom_SignalELS_ChtCoh = sTzoom_OtherNuCC_ChtCoh->ToTH1( pot, colorwheel_mode[5] ); hTzoom_SignalELS_ChtCoh->SetFillColor(colorwheel_mode[5]);
  //
  TH1* hTzoom_Signal_CLECoh    = sTzoom_Signal_CLECoh->ToTH1( pot, colorwheel[0] );         hTzoom_Signal_CLECoh->SetFillColor(colorwheel[0]);
  TH1* hTzoom_OtherNuCC_CLECoh = sTzoom_OtherNuCC_CLECoh->ToTH1( pot, colorwheel[1] );      hTzoom_OtherNuCC_CLECoh->SetFillColor(colorwheel[1]);
  TH1* hTzoom_OtherNuNC_CLECoh = sTzoom_OtherNuNC_CLECoh->ToTH1( pot, colorwheel[2] );      hTzoom_OtherNuNC_CLECoh->SetFillColor(colorwheel[2]);
  TH1* hTzoom_Cosmic_CLECoh    = sTzoom_Cosmic_CLECoh->ToTH1( pot, colorwheel[3] );         hTzoom_Cosmic_CLECoh->SetFillColor(colorwheel[3]);
  TH1* hTzoom_InTime_CLECoh    = sTzoom_InTime_CLECoh->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hTzoom_InTime_CLECoh->SetFillColor(colorwheel[4]);
  TH1* hTzoom_SignalQEL_CLECoh = sTzoom_SignalQEL_CLECoh->ToTH1( pot, colorwheel_mode[0] ); hTzoom_SignalQEL_CLECoh->SetFillColor(colorwheel_mode[0]);
  TH1* hTzoom_SignalMEC_CLECoh = sTzoom_SignalMEC_CLECoh->ToTH1( pot, colorwheel_mode[1] ); hTzoom_SignalMEC_CLECoh->SetFillColor(colorwheel_mode[1]);
  TH1* hTzoom_SignalRES_CLECoh = sTzoom_SignalRES_CLECoh->ToTH1( pot, colorwheel_mode[2] ); hTzoom_SignalRES_CLECoh->SetFillColor(colorwheel_mode[2]);
  TH1* hTzoom_SignalDIS_CLECoh = sTzoom_SignalDIS_CLECoh->ToTH1( pot, colorwheel_mode[3] ); hTzoom_SignalDIS_CLECoh->SetFillColor(colorwheel_mode[3]);
  TH1* hTzoom_SignalCOH_CLECoh = sTzoom_SignalCOH_CLECoh->ToTH1( pot, colorwheel_mode[4] ); hTzoom_SignalCOH_CLECoh->SetFillColor(colorwheel_mode[4]);
  TH1* hTzoom_SignalELS_CLECoh = sTzoom_OtherNuCC_CLECoh->ToTH1( pot, colorwheel_mode[5] ); hTzoom_SignalELS_CLECoh->SetFillColor(colorwheel_mode[5]);

  TH1* hSignal_SelRes    = sSignal_SelRes->ToTH1( pot, colorwheel[0] );         hSignal_SelRes->SetFillColor(colorwheel[0]);
  TH1* hOtherNuCC_SelRes = sOtherNuCC_SelRes->ToTH1( pot, colorwheel[1] );      hOtherNuCC_SelRes->SetFillColor(colorwheel[1]);
  TH1* hOtherNuNC_SelRes = sOtherNuNC_SelRes->ToTH1( pot, colorwheel[2] );      hOtherNuNC_SelRes->SetFillColor(colorwheel[2]);
  TH1* hCosmic_SelRes    = sCosmic_SelRes->ToTH1( pot, colorwheel[3] );         hCosmic_SelRes->SetFillColor(colorwheel[3]);
  TH1* hInTime_SelRes    = sInTime_SelRes->ToTH1( cosmicLivetime, colorwheel[4], kSolid, kLivetime ); hInTime_SelRes->SetFillColor(colorwheel[4]);
  TH1* hSignalQEL_SelRes = sSignalQEL_SelRes->ToTH1( pot, colorwheel_mode[0] ); hSignalQEL_SelRes->SetFillColor(colorwheel_mode[0]);
  TH1* hSignalMEC_SelRes = sSignalMEC_SelRes->ToTH1( pot, colorwheel_mode[1] ); hSignalMEC_SelRes->SetFillColor(colorwheel_mode[1]);
  TH1* hSignalRES_SelRes = sSignalRES_SelRes->ToTH1( pot, colorwheel_mode[2] ); hSignalRES_SelRes->SetFillColor(colorwheel_mode[2]);
  TH1* hSignalDIS_SelRes = sSignalDIS_SelRes->ToTH1( pot, colorwheel_mode[3] ); hSignalDIS_SelRes->SetFillColor(colorwheel_mode[3]);
  TH1* hSignalCOH_SelRes = sSignalCOH_SelRes->ToTH1( pot, colorwheel_mode[4] ); hSignalCOH_SelRes->SetFillColor(colorwheel_mode[4]);
  TH1* hSignalELS_SelRes = sOtherNuCC_SelRes->ToTH1( pot, colorwheel_mode[5] ); hSignalELS_SelRes->SetFillColor(colorwheel_mode[5]);

  // colorwheel_class for the colors here
  TH1* hSignalwProton_recoMuMom_TrueSignal = sSignalwProton_recoMuMom_TrueSignal->ToTH1( pot, colorwheel_class[0] );
  hSignalwProton_recoMuMom_TrueSignal->SetFillColor(colorwheel_class[0]);
  TH1* hSignalwProton_recoMuMom_NumuMuWrong = sSignalwProton_recoMuMom_NumuMuWrong->ToTH1( pot, colorwheel_class[1] );
  hSignalwProton_recoMuMom_NumuMuWrong->SetFillColor(colorwheel_class[1]);
  TH1* hSignalwProton_recoMuMom_NumuPWrong = sSignalwProton_recoMuMom_NumuPWrong->ToTH1( pot, colorwheel_class[2] );
  hSignalwProton_recoMuMom_NumuPWrong->SetFillColor(colorwheel_class[2]);
  TH1* hSignalwProton_recoMuMom_NumuBothWrong = sSignalwProton_recoMuMom_NumuBothWrong->ToTH1( pot, colorwheel_class[3] );
  hSignalwProton_recoMuMom_NumuBothWrong->SetFillColor(colorwheel_class[3]);
  TH1* hSignalwProton_recoMuMom_NumuNonFidVol = sSignalwProton_recoMuMom_NumuNonFidVol->ToTH1( pot, colorwheel_class[4] );
  hSignalwProton_recoMuMom_NumuNonFidVol->SetFillColor(colorwheel_class[4]);
  TH1* hSignalwProton_recoMuMom_NumuOther = sSignalwProton_recoMuMom_NumuOther->ToTH1( pot, colorwheel_class[5] );
  hSignalwProton_recoMuMom_NumuOther->SetFillColor(colorwheel_class[5]);
  TH1* hSignalwProton_recoMuMom_NuOther = sSignalwProton_recoMuMom_NuOther->ToTH1( pot, colorwheel_class[6] );
  hSignalwProton_recoMuMom_NuOther->SetFillColor(colorwheel_class[6]);
  TH1* hSignalwProton_recoMuMom_Cosmic = sSignalwProton_recoMuMom_Cosmic->ToTH1( pot, colorwheel_class[7] );
  hSignalwProton_recoMuMom_Cosmic->SetFillColor(colorwheel_class[7]);
  TH1* hSignalwProton_recoMuMom_Cosmic_InTime = sSignalwProton_recoMuMom_Cosmic_InTime->ToTH1( cosmicLivetime, colorwheel_class[8], kSolid, kLivetime );
  hSignalwProton_recoMuMom_Cosmic_InTime->SetFillColor(colorwheel_class[8]);

  hAll_NoCut->Add(hInTime_NoCut);
  hAll_Selct->Add(hInTime_Selct);

  hAll_P_Selct->Add(hInTime_P_Selct);
  hBackgd_P_Selct->Add(hInTime_P_Selct);

  hAll_PSimple_Selct->Add(hInTime_PSimple_Selct);
  hBackgd_PSimple_Selct->Add(hInTime_PSimple_Selct);

  hAll_Cont_Selct->Add(hInTime_Cont_Selct);
  hAll_SelctG->Add(hInTime_SelctG);


  hAll_P_SelctG->Add(hInTime_P_SelctG);
  hBackgd_P_SelctG->Add(hInTime_P_SelctG);


  hAll_SelTot->Add(hInTime_SelTot);
  hAll_SelCoh->Add(hInTime_SelCoh);
  hT_All_SelCoh->Add(hT_InTime_SelCoh);
  hTzoom_All_SelCoh->Add(hTzoom_InTime_SelCoh);

  hTzoom_All_ChtCoh->Add(hTzoom_InTime_ChtCoh);
  //
  hTzoom_All_CLECoh->Add(hTzoom_InTime_CLECoh);

  hAll_SelRes->Add(hInTime_SelRes);

  hSignalwProton_recoMuMom_BaseSliceSel->Add(hSignalwProton_recoMuMom_Cosmic_InTime);

  double sigEff_NoCut = 1.;
  double sigEff_Selct = hSignal_Selct->Integral() / hSignal_NoCut->Integral();
  double sigEff_Cont_Selct = hSignal_Cont_Selct->Integral() / hSignal_Cont_NoCut->Integral();
  double sigEff_SelctG = hSignal_SelctG->Integral() / hSignal_NoCut->Integral();
  double sigEff_SelTot = hSignal_SelTot->Integral() / hSignal_NoCut->Integral();
  double cohEff_SelCoh = hSignalCOH_SelCoh->Integral() / hSignalCOH_NoCut->Integral();

  // NO CUT
  double totalCounts_NoCut = hAll_NoCut->Integral();
  TLegend *tL_NoCut = new TLegend(0.6,0.6,0.87,0.87);
  tL_NoCut->SetHeader(TString::Format("No Cut = %.1f, Eff: %.1f%%",hSignal_NoCut->Integral(),100.*sigEff_NoCut), "C");
  tL_NoCut->AddEntry(hSignal_NoCut,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_NoCut->AddEntry(hOtherNuCC_NoCut,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_NoCut->AddEntry(hOtherNuNC_NoCut,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_NoCut->AddEntry(hCosmic_NoCut,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_NoCut->AddEntry(hInTime_NoCut,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_NoCut->Integral()/totalCounts_NoCut),"f");

  THStack *hStack_NoCut = new THStack("hStack_NoCut", "");
  hStack_NoCut->Add(hInTime_NoCut);
  hStack_NoCut->Add(hCosmic_NoCut);
  hStack_NoCut->Add(hSignal_NoCut);
  hStack_NoCut->Add(hOtherNuNC_NoCut);
  hStack_NoCut->Add(hOtherNuCC_NoCut);

  new TCanvas;
  hAll_NoCut->Draw("hist");
  hStack_NoCut->Draw("hist same");
  hAll_NoCut->Draw("hist same");
  tL_NoCut->Draw();
  gPad->Print("Spectrum_NoCut.pdf");

  // --> contained version
  TLegend *tL_Cont_NoCut = new TLegend(0.6,0.6,0.87,0.87);
  tL_Cont_NoCut->SetHeader(TString::Format("No Cut = %.1f, Eff: %.1f%%",hSignal_NoCut->Integral(),100.*sigEff_NoCut), "C");
  tL_Cont_NoCut->AddEntry(hSignal_Cont_NoCut,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_Cont_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_Cont_NoCut->AddEntry(hOtherNuCC_Cont_NoCut,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_Cont_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_Cont_NoCut->AddEntry(hOtherNuNC_NoCut,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_Cont_NoCut->AddEntry(hCosmic_NoCut,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_NoCut->Integral()/totalCounts_NoCut),"f");
  tL_Cont_NoCut->AddEntry(hInTime_NoCut,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_NoCut->Integral()/totalCounts_NoCut),"f");

  THStack *hStack_Cont_NoCut = new THStack("hStack_Cont_NoCut", "");
  hStack_Cont_NoCut->Add(hInTime_NoCut);
  hStack_Cont_NoCut->Add(hCosmic_NoCut);
  hStack_Cont_NoCut->Add(hSignal_Cont_NoCut);
  hStack_Cont_NoCut->Add(hOtherNuNC_NoCut);
  hStack_Cont_NoCut->Add(hOtherNuCC_Cont_NoCut);

  new TCanvas;
  hAll_NoCut->Draw("hist");
  hStack_Cont_NoCut->Draw("hist same");
  hAll_NoCut->Draw("hist same");
  tL_Cont_NoCut->Draw();
  gPad->Print("Spectrum_Cont_NoCut.pdf");

  // WITH CUTS
  double totalCounts_Selct = hAll_Selct->Integral();
  TLegend *tL_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_Selct->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_Selct->Integral(),100.*sigEff_Selct), "C");
  tL_Selct->AddEntry(hSignal_Selct,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_Selct->Integral()/totalCounts_Selct),"f");
  tL_Selct->AddEntry(hOtherNuCC_Selct,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_Selct->Integral()/totalCounts_Selct),"f");
  tL_Selct->AddEntry(hOtherNuNC_Selct,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_Selct->Integral()/totalCounts_Selct),"f");
  tL_Selct->AddEntry(hCosmic_Selct,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_Selct->Integral()/totalCounts_Selct),"f");
  tL_Selct->AddEntry(hInTime_Selct,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_Selct->Integral()/totalCounts_Selct),"f");

  THStack *hStack_Selct = new THStack("hStack_Selct", "");
  hStack_Selct->Add(hSignal_Selct);
  hStack_Selct->Add(hOtherNuNC_Selct);
  hStack_Selct->Add(hOtherNuCC_Selct);
  hStack_Selct->Add(hInTime_Selct);
  hStack_Selct->Add(hCosmic_Selct);

  new TCanvas;
  hAll_Selct->Draw("hist");
  hStack_Selct->Draw("hist same");
  hAll_Selct->Draw("hist same");
  tL_Selct->Draw();
  gPad->Print("Spectrum_Selected.pdf");

  // --> that but with modes
  TLegend *tLMode_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_Selct->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_Selct->Integral(),100.*sigEff_Selct), "C");
  tLMode_Selct->AddEntry(hSignalQEL_Selct,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hSignalMEC_Selct,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hSignalRES_Selct,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hSignalDIS_Selct,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hSignalCOH_Selct,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hSignalELS_Selct,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hOtherNuNC_Selct,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hCosmic_Selct,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_Selct->Integral()/totalCounts_Selct),"f");
  tLMode_Selct->AddEntry(hInTime_Selct,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_Selct->Integral()/totalCounts_Selct),"f");

  THStack *hStackMode_Selct = new THStack("hStackMode_Selct", "");
  hStackMode_Selct->Add(hSignalQEL_Selct);
  hStackMode_Selct->Add(hSignalMEC_Selct);
  hStackMode_Selct->Add(hSignalRES_Selct);
  hStackMode_Selct->Add(hSignalDIS_Selct);
  hStackMode_Selct->Add(hSignalCOH_Selct);
  hStackMode_Selct->Add(hSignalELS_Selct);
  hStackMode_Selct->Add(hOtherNuNC_Selct);
  hStackMode_Selct->Add(hInTime_Selct);
  hStackMode_Selct->Add(hCosmic_Selct);

  new TCanvas;
  hAll_Selct->Draw("hist");
  hStackMode_Selct->Draw("hist same");
  hAll_Selct->Draw("hist same");
  tLMode_Selct->Draw();
  gPad->Print("Spectrum_Modes_Selected.pdf");

  // PROTONS
  double totalCounts_P_Selct = hAll_P_Selct->Integral();
  TLegend *tL_P_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_P_Selct->SetHeader(TString::Format("Select slices = %.1f",hSignal_P_Selct->Integral()), "C");
  tL_P_Selct->AddEntry(hSignal_P_Selct,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_P_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_P_Selct->AddEntry(hOtherNuCC_P_Selct,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_P_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_P_Selct->AddEntry(hOtherNuNC_P_Selct,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_P_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_P_Selct->AddEntry(hCosmic_P_Selct,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_P_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_P_Selct->AddEntry(hInTime_P_Selct,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_P_Selct->Integral()/totalCounts_P_Selct),"f");

  THStack *hStack_P_Selct = new THStack("hStack_P_Selct", "");
  hStack_P_Selct->Add(hSignal_P_Selct);
  hStack_P_Selct->Add(hOtherNuNC_P_Selct);
  hStack_P_Selct->Add(hOtherNuCC_P_Selct);
  hStack_P_Selct->Add(hInTime_P_Selct);
  hStack_P_Selct->Add(hCosmic_P_Selct);

  new TCanvas;
  hAll_P_Selct->Draw("hist");
  hStack_P_Selct->Draw("hist same");
  hAll_P_Selct->Draw("hist same");
  tL_P_Selct->Draw();
  gPad->Print("Spectrum_Selected_LeadingPCand.pdf");

  // true p vs not true p vs bckgd
  TLegend *tL_TruP_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_TruP_Selct->SetHeader(TString::Format("Select slices = %.1f",hSignal_P_Selct->Integral()), "C");
  tL_TruP_Selct->AddEntry(hSignal_TruP_Selct,TString::Format("Signal True p: %.1f%%",100.*hSignal_TruP_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_TruP_Selct->AddEntry(hSignal_NotP_Selct,TString::Format("Signal MisID: %.1f%%", 100.*hSignal_NotP_Selct->Integral()/totalCounts_P_Selct),"f");
  tL_TruP_Selct->AddEntry(hBackgd_P_Selct,   TString::Format("Backgrounds: %.1f%%",  100.*hBackgd_P_Selct->Integral()/totalCounts_P_Selct),   "f");

  THStack *hStack_TruP_Selct = new THStack("hStack_TruP_Selct", "");
  hStack_TruP_Selct->Add(hSignal_TruP_Selct);
  hStack_TruP_Selct->Add(hSignal_NotP_Selct);
  hStack_TruP_Selct->Add(hBackgd_P_Selct);

  new TCanvas;
  hAll_P_Selct->Draw("hist");
  hStack_TruP_Selct->Draw("hist same");
  hAll_P_Selct->Draw("hist same");
  tL_TruP_Selct->Draw();
  gPad->Print("Spectrum_Selected_LeadingPCand_TruP.pdf");

  // p momentum reco vs true
  new TCanvas;
  hRecoVsTrue_ProtonP->SetTitle("");
  hRecoVsTrue_ProtonP->Draw("colz");
  TLine *tline = new TLine(0.0,0.0,1.6,1.6);
  tline->SetLineWidth(2);
  tline->SetLineColor(kGray+1);
  tline->SetLineStyle(2);
  tline->Draw("same");
  gPad->Print("Spectrum_Selected_LeadingPCand_TruP_RecoVsTrue.pdf");


  // PROTONS with Selection minus the proton cut...
  double totalCounts_PSimple_Selct = hAll_PSimple_Selct->Integral();
  TLegend *tL_PSimple_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_PSimple_Selct->SetHeader(TString::Format("Select slices = %.1f",hSignal_PSimple_Selct->Integral()), "C");
  tL_PSimple_Selct->AddEntry(hSignal_PSimple_Selct,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_PSimple_Selct->AddEntry(hOtherNuCC_PSimple_Selct,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_PSimple_Selct->AddEntry(hOtherNuNC_PSimple_Selct,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_PSimple_Selct->AddEntry(hCosmic_PSimple_Selct,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_PSimple_Selct->AddEntry(hInTime_PSimple_Selct,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");

  THStack *hStack_PSimple_Selct = new THStack("hStack_PSimple_Selct", "");
  hStack_PSimple_Selct->Add(hSignal_PSimple_Selct);
  hStack_PSimple_Selct->Add(hOtherNuNC_PSimple_Selct);
  hStack_PSimple_Selct->Add(hOtherNuCC_PSimple_Selct);
  hStack_PSimple_Selct->Add(hInTime_PSimple_Selct);
  hStack_PSimple_Selct->Add(hCosmic_PSimple_Selct);

  new TCanvas;
  hAll_PSimple_Selct->Draw("hist");
  hStack_PSimple_Selct->Draw("hist same");
  hAll_PSimple_Selct->Draw("hist same");
  tL_PSimple_Selct->Draw();
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE.pdf");

  // -- true p vs not true p vs bckgd
  TLegend *tL_TruPSimple_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_TruPSimple_Selct->SetHeader(TString::Format("Select slices = %.1f",hSignal_PSimple_Selct->Integral()), "C");
  tL_TruPSimple_Selct->AddEntry(hSignal_TruPSimple_Selct,TString::Format("Signal True p: %.1f%%",100.*hSignal_TruPSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_TruPSimple_Selct->AddEntry(hSignal_NotPSimple_Selct,TString::Format("Signal MisID: %.1f%%", 100.*hSignal_NotPSimple_Selct->Integral()/totalCounts_PSimple_Selct),"f");
  tL_TruPSimple_Selct->AddEntry(hBackgd_PSimple_Selct,   TString::Format("Backgrounds: %.1f%%",  100.*hBackgd_PSimple_Selct->Integral()/totalCounts_PSimple_Selct),   "f");

  THStack *hStack_TruPSimple_Selct = new THStack("hStack_TruPSimple_Selct", "");
  hStack_TruPSimple_Selct->Add(hSignal_TruPSimple_Selct);
  hStack_TruPSimple_Selct->Add(hSignal_NotPSimple_Selct);
  hStack_TruPSimple_Selct->Add(hBackgd_PSimple_Selct);

  new TCanvas;
  hAll_PSimple_Selct->Draw("hist");
  hStack_TruPSimple_Selct->Draw("hist same");
  hAll_PSimple_Selct->Draw("hist same");
  tL_TruPSimple_Selct->Draw();
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE_TruP.pdf");

  // -- chi2 proton
  new TCanvas;
  hSignal_TruPSimple_Selct_Chi2P->Draw("hist");
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE_TruP_Chi2P.pdf");

  new TCanvas;
  hSignal_TruStopPSimple_Selct_Chi2P->Draw("hist");
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE_Stopping_TruP_Chi2P.pdf");

  // -- chi2 muon
  new TCanvas;
  hSignal_TruPSimple_Selct_Chi2Mu->Draw("hist");
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE_TruP_Chi2Mu.pdf");

  new TCanvas;
  hSignal_TruStopPSimple_Selct_Chi2Mu->Draw("hist");
  gPad->Print("Spectrum_Selected_LeadingPCand_SIMPLE_Stopping_TruP_Chi2Mu.pdf");



  // proton multiplicities
  // Plot 1 no cut
  new TCanvas;
  TLegend *tL_PMultNoCut = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultNoCut->AddEntry( hSignal_PMultPreFSI_NoCut, "Pre-FSI", "l" );
  tL_PMultNoCut->AddEntry( hSignal_PMultPostFSI_NoCut, "Post-FSI, No Cut", "l" );
  tL_PMultNoCut->AddEntry( hSignal_PMultPostFSIMinE_NoCut, "Post, E>15MeV", "l" );
  tL_PMultNoCut->AddEntry( hSignal_PMultPostFSIE100_NoCut, "Post, E>100MeV", "l" );

  hSignal_PMultPreFSI_NoCut->Draw("hist");
  hSignal_PMultPostFSI_NoCut->Draw("hist same");
  hSignal_PMultPostFSIMinE_NoCut->Draw("hist same");
  hSignal_PMultPostFSIE100_NoCut->Draw("hist same");
  tL_PMultNoCut->Draw();
  gPad->Print("Spectrum_NoCut_ProtonMultiplicity_True.pdf");
  
  // Plot 1 w/ cut
  new TCanvas;
  TLegend *tL_PMultSelct = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultSelct->AddEntry( hSignal_PMultPreFSI_Selct, "Pre-FSI", "l" );
  tL_PMultSelct->AddEntry( hSignal_PMultPostFSI_Selct, "Post-FSI, No Cut", "l" );
  tL_PMultSelct->AddEntry( hSignal_PMultPostFSIMinE_Selct, "Post, E>15MeV", "l" );
  tL_PMultSelct->AddEntry( hSignal_PMultPostFSIE100_Selct, "Post, E>100MeV", "l" );

  hSignal_PMultPreFSI_Selct->Draw("hist");
  hSignal_PMultPostFSI_Selct->Draw("hist same");
  hSignal_PMultPostFSIMinE_Selct->Draw("hist same");
  hSignal_PMultPostFSIE100_Selct->Draw("hist same");
  tL_PMultSelct->Draw();
  gPad->Print("Spectrum_Selct_ProtonMultiplicity_True.pdf");

  // Plot 2 no cut
  new TCanvas;
  TLegend *tL_PMultNoCutReco = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultNoCutReco->AddEntry( hSignal_PMultPostFSIE100_RecoComp_NoCut, "Post-FSI, E>100MeV", "l" );
  tL_PMultNoCutReco->AddEntry( hSignal_PMultReco_NoCut, "Reco w/o stubs", "l" );
  tL_PMultNoCutReco->AddEntry( hSignal_PMultRecoCheated_NoCut, "'Cheated'", "l" );

  hSignal_PMultReco_NoCut->Draw("hist");
  hSignal_PMultPostFSIE100_RecoComp_NoCut->Draw("hist same");
  hSignal_PMultRecoCheated_NoCut->Draw("hist same");
  tL_PMultNoCutReco->Draw();
  gPad->Print("Spectrum_NoCut_ProtonMultiplicity_Reco.pdf");
  
  // Plot 2 w/ cut
  new TCanvas;
  TLegend *tL_PMultSelctReco = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultSelctReco->AddEntry( hSignal_PMultPostFSIE100_RecoComp_Selct, "Post-FSI, E>100MeV", "l" );
  tL_PMultSelctReco->AddEntry( hSignal_PMultReco_Selct, "Reco w/o stubs", "l" );
  tL_PMultSelctReco->AddEntry( hSignal_PMultRecoCheated_Selct, "'Cheated'", "l" );

  hSignal_PMultReco_Selct->Draw("hist");
  hSignal_PMultPostFSIE100_RecoComp_Selct->Draw("hist same");
  hSignal_PMultRecoCheated_Selct->Draw("hist same");
  tL_PMultSelctReco->Draw();
  gPad->Print("Spectrum_Selct_ProtonMultiplicity_Reco.pdf");

  // Plot 3 no cut
  new TCanvas;
  TLegend *tL_PMultNoCutRecoStubs = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultNoCutRecoStubs->AddEntry( hSignal_PMultPostFSIMinE_RecoComp_NoCut, "Post-FSI, E>15MeV", "l" );
  tL_PMultNoCutRecoStubs->AddEntry( hSignal_PMultRecoStubs_NoCut, "Reco w/ stubs", "l" );

  hSignal_PMultPostFSIMinE_RecoComp_NoCut->Draw("hist");
  hSignal_PMultRecoStubs_NoCut->Draw("hist same");
  tL_PMultNoCutRecoStubs->Draw();
  gPad->Print("Spectrum_NoCut_ProtonMultiplicity_RecoStubs.pdf");

  // Plot 3 w/ cut
  new TCanvas;
  TLegend *tL_PMultSelctRecoStubs = new TLegend(0.6,0.6,0.87,0.87);
  tL_PMultSelctRecoStubs->AddEntry( hSignal_PMultPostFSIMinE_RecoComp_Selct, "Post-FSI, E>15MeV", "l" );
  tL_PMultSelctRecoStubs->AddEntry( hSignal_PMultRecoStubs_Selct, "Reco w/ stubs", "l" );

  hSignal_PMultPostFSIMinE_RecoComp_Selct->Draw("hist");
  hSignal_PMultRecoStubs_Selct->Draw("hist same");
  tL_PMultSelctRecoStubs->Draw();
  gPad->Print("Spectrum_Selct_ProtonMultiplicity_RecoStubs.pdf");

  // True multiplicity/cheated efficiency...
  new TCanvas;
  TLegend *tL_TrueProtonKE = new TLegend(0.6,0.6,0.87,0.87);
  tL_TrueProtonKE->AddEntry( hSignal_TrueProtonKE_SpillMV, "True KE distribution", "l" );
  tL_TrueProtonKE->AddEntry( hSignal_TrueProtonKE_CheatedReco_SpillMV, "True KE w reco match", "l" );

  hSignal_TrueProtonKE_SpillMV->Draw("hist");
  hSignal_TrueProtonKE_CheatedReco_SpillMV->Draw("hist same");
  tL_TrueProtonKE->Draw();
  gPad->Print("Spectrum_TrueProtonKE.pdf");

  new TCanvas;
  rSignal_TrueProtonKE_SpillMV->Draw("");
  gPad->Print("Spectrum_TrueProtonKE_CheatedRecoEff.pdf");

  new TCanvas;
  TLegend *tL_TrueProtonKE_RES = new TLegend(0.6,0.6,0.87,0.87);
  tL_TrueProtonKE_RES->AddEntry( hSignal_TrueProtonKE_RES_SpillMV, "True KE distribution", "l" );
  tL_TrueProtonKE_RES->AddEntry( hSignal_TrueProtonKE_RES_CheatedReco_SpillMV, "True KE w reco match", "l" );

  hSignal_TrueProtonKE_RES_SpillMV->Draw("hist");
  hSignal_TrueProtonKE_RES_CheatedReco_SpillMV->Draw("hist same");
  tL_TrueProtonKE_RES->Draw();
  gPad->Print("Spectrum_TrueProtonKE_RES.pdf");

  new TCanvas;
  rSignal_TrueProtonKE_RES_SpillMV->Draw("");
  gPad->Print("Spectrum_TrueProtonKE_RES_CheatedRecoEff.pdf");

  new TCanvas;
  TLegend *tL_TrueProtonMom = new TLegend(0.6,0.6,0.87,0.87);
  tL_TrueProtonMom->AddEntry( hSignal_TrueProtonMom_SpillMV, "True |p| distribution", "l" );
  tL_TrueProtonMom->AddEntry( hSignal_TrueProtonMom_CheatedReco_SpillMV, "True |p| w reco match", "l" );

  hSignal_TrueProtonMom_SpillMV->Draw("hist");
  hSignal_TrueProtonMom_CheatedReco_SpillMV->Draw("hist same");
  tL_TrueProtonMom->Draw();
  gPad->Print("Spectrum_TrueProtonMom.pdf");

  new TCanvas;
  rSignal_TrueProtonMom_SpillMV->Draw("");
  gPad->Print("Spectrum_TrueProtonMom_CheatedRecoEff.pdf");

  // neutrino E
  new TCanvas;
  hSignal_TrueNuE_SpillMV->Draw("hist");
  hSignal_TrueNuE_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_Signal_TrueNuE_Selct.pdf");

  new TCanvas;
  rSignal_TrueNuE_SpillMV->Draw("");
  gPad->Print("Spectrum_Signal_TrueNuE_Selct_Eff.pdf");

  new TCanvas;
  hSignalwProton_TrueMuMom_SpillMV->Draw("hist");
  hSignalwProton_TrueMuMom_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_SignalwProton_TrueMuMom_Selct.pdf");

  new TCanvas;
  hAll_RecoMuMom_Selct_SpillMV->Draw("hist");
  hSignalwProton_RecoMuMom_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_SignalAndNotSignal_RecoMuMom_Selct.pdf");

  new TCanvas;
  rSignalwProton_TrueMuMom_SpillMV->Draw("");
  gPad->Print("Spectrum_SignalwProton_TrueMuMom_Selct_Eff.pdf");
  
  new TCanvas;
  rSignalwProton_RecoMuMom_Selct_SpillMV->Draw("");
  gPad->Print("Spectrum_SignalwProton_RecoMuMom_Selct_Pur.pdf");

  // Angles...
  new TCanvas;
  hNuMIMuons_CosThXZ->Draw("hist");
  gPad->Print("Spectrum_MuFromNu_CosThXZ.pdf");

  new TCanvas;
  hNuMIMuons_CosThNuMI->Draw("hist");
  hNuMIMuons_CosThNuMINotIsoch->Draw("hist same");
  hNuMIMuons_CosThNuMINotIsochNotPerp->Draw("hist same");
  gPad->Print("Spectrum_MuFromNu_CosThNuMI.pdf");

  // L and L/E
  new TCanvas;
  hSignal_TrueNuL_SpillMV->Draw("hist");
  //hSignal_TrueNuL_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_Signal_TrueNuL_Selct.pdf");

  new TCanvas;
  hSignal_TrueNuLOverE_SpillMV->Draw("hist");
  //hSignal_TrueNuLOverE_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_Signal_TrueNuLOverE_Selct.pdf");

  new TCanvas;
  hNumuCC_TrueNuE_SpillMV->Draw("hist");
  //hNumuCC_TrueNuE_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_NumuCC_TrueNuE_Selct.pdf");

  new TCanvas;
  hNumuCC_TrueNuL_SpillMV->Draw("hist");
  //hNumuCC_TrueNuL_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_NumuCC_TrueNuL_Selct.pdf");

  new TCanvas;
  hNumuCC_TrueNuLOverE_SpillMV->Draw("hist");
  //hNumuCC_TrueNuLOverE_Selct_SpillMV->Draw("hist same");
  gPad->Print("Spectrum_NumuCC_TrueNuLOverE_Selct.pdf");


  
  // SOME 2D VERSIONS OF MULTIPLICITY
  // 2d: sRecoVsTrue_ProtonMult_NoCut
  new TCanvas;
  hRecoVsTrue_ProtonMult_NoCut->SetTitle("");
  hRecoVsTrue_ProtonMult_NoCut->Draw("colz");
  TLine *tline2 = new TLine(-0.5,-0.5,6.5,6.5);
  tline2->SetLineWidth(2);
  tline2->SetLineColor(kGray+1);
  tline2->SetLineStyle(2);
  tline2->Draw("same");
  gPad->Print("Spectrum_NoCut_ProtonMultiplicity2D.pdf");

  new TCanvas;
  hRecoVsTrue_ProtonMult_Selct->SetTitle("");
  hRecoVsTrue_ProtonMult_Selct->Draw("colz");
  TLine *tline3 = new TLine(-0.5,-0.5,6.5,6.5);
  tline3->SetLineWidth(2);
  tline3->SetLineColor(kGray+1);
  tline3->SetLineStyle(2);
  tline3->Draw("same");
  gPad->Print("Spectrum_Selct_ProtonMultiplicity2D.pdf");

  // CONTAINED SELECTION
  double totalCounts_Cont_Selct = hAll_Cont_Selct->Integral();
  TLegend *tL_Cont_Selct = new TLegend(0.6,0.6,0.87,0.87);
  tL_Cont_Selct->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_Cont_Selct->Integral(),100.*sigEff_Cont_Selct), "C");
  tL_Cont_Selct->AddEntry(hSignal_Cont_Selct,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_Cont_Selct->Integral()/totalCounts_Cont_Selct),"f");
  tL_Cont_Selct->AddEntry(hOtherNuCC_Cont_Selct,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_Cont_Selct->Integral()/totalCounts_Cont_Selct),"f");
  tL_Cont_Selct->AddEntry(hOtherNuNC_Cont_Selct,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_Cont_Selct->Integral()/totalCounts_Cont_Selct),"f");
  tL_Cont_Selct->AddEntry(hCosmic_Cont_Selct,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_Cont_Selct->Integral()/totalCounts_Cont_Selct),"f");
  tL_Cont_Selct->AddEntry(hInTime_Cont_Selct,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_Cont_Selct->Integral()/totalCounts_Cont_Selct),"f");

  THStack *hStack_Cont_Selct = new THStack("hStack_Cont_Selct", "");
  hStack_Cont_Selct->Add(hSignal_Cont_Selct);
  hStack_Cont_Selct->Add(hOtherNuNC_Cont_Selct);
  hStack_Cont_Selct->Add(hOtherNuCC_Cont_Selct);
  hStack_Cont_Selct->Add(hInTime_Cont_Selct);
  hStack_Cont_Selct->Add(hCosmic_Cont_Selct);

  new TCanvas;
  hAll_Cont_Selct->Draw("hist");
  hStack_Cont_Selct->Draw("hist same");
  hAll_Cont_Selct->Draw("hist same");
  tL_Cont_Selct->Draw();
  gPad->Print("Spectrum_Cont_Selected.pdf");


  // Gray's proposed sample
  double totalCounts_SelctG = hAll_SelctG->Integral();
  TLegend *tL_SelctG = new TLegend(0.6,0.6,0.87,0.87);
  tL_SelctG->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_SelctG->Integral(),100.*sigEff_SelctG), "C");
  tL_SelctG->AddEntry(hSignal_SelctG,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_SelctG->Integral()/totalCounts_SelctG),"f");
  tL_SelctG->AddEntry(hOtherNuCC_SelctG,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_SelctG->Integral()/totalCounts_SelctG),"f");
  tL_SelctG->AddEntry(hOtherNuNC_SelctG,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelctG->Integral()/totalCounts_SelctG),"f");
  tL_SelctG->AddEntry(hCosmic_SelctG,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelctG->Integral()/totalCounts_SelctG),"f");
  tL_SelctG->AddEntry(hInTime_SelctG,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelctG->Integral()/totalCounts_SelctG),"f");

  THStack *hStack_SelctG = new THStack("hStack_SelctG", "");
  hStack_SelctG->Add(hSignal_SelctG);
  hStack_SelctG->Add(hOtherNuNC_SelctG);
  hStack_SelctG->Add(hOtherNuCC_SelctG);
  hStack_SelctG->Add(hInTime_SelctG);
  hStack_SelctG->Add(hCosmic_SelctG);

  new TCanvas;
  hAll_SelctG->Draw("hist");
  hStack_SelctG->Draw("hist same");
  hAll_SelctG->Draw("hist same");
  tL_SelctG->Draw();
  gPad->Print("Spectrum_GraysProposedSample.pdf");

  // Version with interaction modes
  TLegend *tLMode_SelctG = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_SelctG->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_SelctG->Integral(),100.*sigEff_SelctG), "C");
  tLMode_SelctG->AddEntry(hSignalQEL_SelctG,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hSignalMEC_SelctG,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hSignalRES_SelctG,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hSignalDIS_SelctG,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hSignalCOH_SelctG,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hSignalELS_SelctG,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hOtherNuNC_SelctG,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hCosmic_SelctG,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelctG->Integral()/totalCounts_SelctG),"f");
  tLMode_SelctG->AddEntry(hInTime_SelctG,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelctG->Integral()/totalCounts_SelctG),"f");

  THStack *hStackMode_SelctG = new THStack("hStackMode_SelctG", "");
  hStackMode_SelctG->Add(hSignalQEL_SelctG);
  hStackMode_SelctG->Add(hSignalMEC_SelctG);
  hStackMode_SelctG->Add(hSignalRES_SelctG);
  hStackMode_SelctG->Add(hSignalDIS_SelctG);
  hStackMode_SelctG->Add(hSignalCOH_SelctG);
  hStackMode_SelctG->Add(hSignalELS_SelctG);
  hStackMode_SelctG->Add(hOtherNuNC_SelctG);
  hStackMode_SelctG->Add(hInTime_SelctG);
  hStackMode_SelctG->Add(hCosmic_SelctG);

  new TCanvas;
  hAll_SelctG->Draw("hist");
  hStackMode_SelctG->Draw("hist same");
  hAll_SelctG->Draw("hist same");
  tLMode_SelctG->Draw();
  gPad->Print("Spectrum_Modes_GraysProposedSample.pdf");


  // PROTONS in Gray's sample
  double totalCounts_P_SelctG = hAll_P_SelctG->Integral();
  TLegend *tL_P_SelctG = new TLegend(0.6,0.6,0.87,0.87);
  tL_P_SelctG->SetHeader(TString::Format("Select slices = %.1f",hSignal_P_SelctG->Integral()), "C");
  tL_P_SelctG->AddEntry(hSignal_P_SelctG,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_P_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_P_SelctG->AddEntry(hOtherNuCC_P_SelctG,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_P_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_P_SelctG->AddEntry(hOtherNuNC_P_SelctG,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_P_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_P_SelctG->AddEntry(hCosmic_P_SelctG,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_P_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_P_SelctG->AddEntry(hInTime_P_SelctG,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_P_SelctG->Integral()/totalCounts_P_SelctG),"f");

  THStack *hStack_P_SelctG = new THStack("hStack_P_SelctG", "");
  hStack_P_SelctG->Add(hSignal_P_SelctG);
  hStack_P_SelctG->Add(hOtherNuNC_P_SelctG);
  hStack_P_SelctG->Add(hOtherNuCC_P_SelctG);
  hStack_P_SelctG->Add(hInTime_P_SelctG);
  hStack_P_SelctG->Add(hCosmic_P_SelctG);

  new TCanvas;
  hAll_P_SelctG->Draw("hist");
  hStack_P_SelctG->Draw("hist same");
  hAll_P_SelctG->Draw("hist same");
  tL_P_SelctG->Draw();
  gPad->Print("Spectrum_GraysSample_LeadingPCand.pdf");

  // true p vs not true p vs bckgd
  TLegend *tL_TruP_SelctG = new TLegend(0.6,0.6,0.87,0.87);
  tL_TruP_SelctG->SetHeader(TString::Format("Select slices = %.1f",hSignal_P_SelctG->Integral()), "C");
  tL_TruP_SelctG->AddEntry(hSignal_TruP_SelctG,TString::Format("Signal True p: %.1f%%",100.*hSignal_TruP_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_TruP_SelctG->AddEntry(hSignal_NotP_SelctG,TString::Format("Signal MisID: %.1f%%", 100.*hSignal_NotP_SelctG->Integral()/totalCounts_P_SelctG),"f");
  tL_TruP_SelctG->AddEntry(hBackgd_P_SelctG,   TString::Format("Backgrounds: %.1f%%",  100.*hBackgd_P_SelctG->Integral()/totalCounts_P_SelctG),   "f");

  THStack *hStack_TruP_SelctG = new THStack("hStack_TruP_SelctG", "");
  hStack_TruP_SelctG->Add(hSignal_TruP_SelctG);
  hStack_TruP_SelctG->Add(hSignal_NotP_SelctG);
  hStack_TruP_SelctG->Add(hBackgd_P_SelctG);

  new TCanvas;
  hAll_P_SelctG->Draw("hist");
  hStack_TruP_SelctG->Draw("hist same");
  hAll_P_SelctG->Draw("hist same");
  tL_TruP_SelctG->Draw();
  gPad->Print("Spectrum_GraysSample_LeadingPCand_TruP.pdf");


  // Gray's proposed sample + Full selection cuts
  double totalCounts_SelTot = hAll_SelTot->Integral();
  TLegend *tL_SelTot = new TLegend(0.6,0.6,0.87,0.87);
  tL_SelTot->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_SelTot->Integral(),100.*sigEff_SelTot), "C");
  tL_SelTot->AddEntry(hSignal_SelTot,TString::Format("NuMu CC Signal: %.1f%%",100.*hSignal_SelTot->Integral()/totalCounts_SelTot),"f");
  tL_SelTot->AddEntry(hOtherNuCC_SelTot,TString::Format("Other Nu CC: %.1f%%",100.*hOtherNuCC_SelTot->Integral()/totalCounts_SelTot),"f");
  tL_SelTot->AddEntry(hOtherNuNC_SelTot,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelTot->Integral()/totalCounts_SelTot),"f");
  tL_SelTot->AddEntry(hCosmic_SelTot,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelTot->Integral()/totalCounts_SelTot),"f");
  tL_SelTot->AddEntry(hInTime_SelTot,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelTot->Integral()/totalCounts_SelTot),"f");

  THStack *hStack_SelTot = new THStack("hStack_SelTot", "");
  hStack_SelTot->Add(hSignal_SelTot);
  hStack_SelTot->Add(hOtherNuNC_SelTot);
  hStack_SelTot->Add(hOtherNuCC_SelTot);
  hStack_SelTot->Add(hInTime_SelTot);
  hStack_SelTot->Add(hCosmic_SelTot);

  new TCanvas;
  hAll_SelTot->Draw("hist");
  hStack_SelTot->Draw("hist same");
  hAll_SelTot->Draw("hist same");
  tL_SelTot->Draw();
  gPad->Print("Spectrum_GraysProposedSample_AndFullSelection.pdf");

  // -- and based on int modes
  TLegend *tLMode_SelTot = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_SelTot->SetHeader(TString::Format("Select = %.1f, Eff: %.1f%%",hSignal_SelTot->Integral(),100.*sigEff_SelTot), "C");
  tLMode_SelTot->AddEntry(hSignalQEL_SelTot,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hSignalMEC_SelTot,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hSignalRES_SelTot,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hSignalDIS_SelTot,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hSignalELS_SelTot,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hSignalCOH_SelTot,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hOtherNuNC_SelTot,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hCosmic_SelTot,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelTot->Integral()/totalCounts_SelTot),"f");
  tLMode_SelTot->AddEntry(hInTime_SelTot,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelTot->Integral()/totalCounts_SelTot),"f");

  THStack *hStackMode_SelTot = new THStack("hStackMode_SelTot", "");
  hStackMode_SelTot->Add(hSignalQEL_SelTot);
  hStackMode_SelTot->Add(hSignalMEC_SelTot);
  hStackMode_SelTot->Add(hSignalRES_SelTot);
  hStackMode_SelTot->Add(hSignalDIS_SelTot);
  hStackMode_SelTot->Add(hSignalCOH_SelTot);
  hStackMode_SelTot->Add(hSignalELS_SelTot);
  hStackMode_SelTot->Add(hOtherNuNC_SelTot);
  hStackMode_SelTot->Add(hInTime_SelTot);
  hStackMode_SelTot->Add(hCosmic_SelTot);

  new TCanvas;
  hAll_SelTot->Draw("hist");
  hStackMode_SelTot->Draw("hist same");
  hAll_SelTot->Draw("hist same");
  tLMode_SelTot->Draw();
  gPad->Print("Spectrum_Modes_GraysProposedSample_AndFullSelection.pdf");


  ////////////////////
  // Coherent study //
  ////////////////////

  // Interaction modes no cuts
  TLegend *tLMode_NoCut = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_NoCut->SetHeader(TString::Format("COH Select = %.1f, Eff: %.1f%%",hSignalCOH_NoCut->Integral(),100.*sigEff_NoCut), "C");
  tLMode_NoCut->AddEntry(hSignalQEL_NoCut,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hSignalMEC_NoCut,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hSignalRES_NoCut,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hSignalDIS_NoCut,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hSignalCOH_NoCut,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hSignalELS_NoCut,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hOtherNuNC_NoCut,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hCosmic_NoCut,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_NoCut->Integral()/totalCounts_NoCut),"f");
  tLMode_NoCut->AddEntry(hInTime_NoCut,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_NoCut->Integral()/totalCounts_NoCut),"f");

  THStack *hStackMode_NoCut = new THStack("hStackMode_NoCut", "");
  hStackMode_NoCut->Add(hSignalCOH_NoCut);
  hStackMode_NoCut->Add(hSignalQEL_NoCut);
  hStackMode_NoCut->Add(hSignalMEC_NoCut);
  hStackMode_NoCut->Add(hSignalRES_NoCut);
  hStackMode_NoCut->Add(hSignalDIS_NoCut);
  hStackMode_NoCut->Add(hSignalELS_NoCut);
  hStackMode_NoCut->Add(hOtherNuNC_NoCut);
  hStackMode_NoCut->Add(hInTime_NoCut);
  hStackMode_NoCut->Add(hCosmic_NoCut);

  new TCanvas;
  hAll_NoCut->Draw("hist");
  hStackMode_NoCut->Draw("hist same");
  hAll_NoCut->Draw("hist same");
  tLMode_NoCut->Draw();
  gPad->Print("Spectrum_Modes_CoherentStudy_NoCut.pdf");

  // With cuts
  double totalCounts_SelCoh = hAll_SelCoh->Integral();
  TLegend *tLMode_SelCoh = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_SelCoh->SetHeader(TString::Format("COH Select = %.1f, Eff: %.1f%%",hSignalCOH_SelCoh->Integral(),100.*cohEff_SelCoh), "C");
  tLMode_SelCoh->AddEntry(hSignalQEL_SelCoh,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hSignalMEC_SelCoh,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hSignalRES_SelCoh,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hSignalDIS_SelCoh,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hSignalCOH_SelCoh,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hSignalELS_SelCoh,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hOtherNuNC_SelCoh,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hCosmic_SelCoh,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelCoh->Integral()/totalCounts_SelCoh),"f");
  tLMode_SelCoh->AddEntry(hInTime_SelCoh,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelCoh->Integral()/totalCounts_SelCoh),"f");

  THStack *hStackMode_SelCoh = new THStack("hStackMode_SelCoh", "");
  hStackMode_SelCoh->Add(hSignalCOH_SelCoh);
  hStackMode_SelCoh->Add(hSignalQEL_SelCoh);
  hStackMode_SelCoh->Add(hSignalMEC_SelCoh);
  hStackMode_SelCoh->Add(hSignalRES_SelCoh);
  hStackMode_SelCoh->Add(hSignalDIS_SelCoh);
  hStackMode_SelCoh->Add(hSignalELS_SelCoh);
  hStackMode_SelCoh->Add(hOtherNuNC_SelCoh);
  hStackMode_SelCoh->Add(hInTime_SelCoh);
  hStackMode_SelCoh->Add(hCosmic_SelCoh);

  new TCanvas;
  hAll_SelCoh->Draw("hist");
  hStackMode_SelCoh->Draw("hist same");
  hAll_SelCoh->Draw("hist same");
  tLMode_SelCoh->Draw();
  gPad->Print("Spectrum_Modes_CoherentStudy_SelCoh.pdf");

  // |t|
  TLegend *tLMode_T_SelCoh = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_T_SelCoh->SetHeader("|t| in COH Selection","C");
  tLMode_T_SelCoh->AddEntry(hT_SignalQEL_SelCoh,"NuMu CC QEL","f");
  tLMode_T_SelCoh->AddEntry(hT_SignalMEC_SelCoh,"NuMu CC MEC","f");
  tLMode_T_SelCoh->AddEntry(hT_SignalRES_SelCoh,"NuMu CC RES","f");
  tLMode_T_SelCoh->AddEntry(hT_SignalDIS_SelCoh,"NuMu CC DIS","f");
  tLMode_T_SelCoh->AddEntry(hT_SignalCOH_SelCoh,"NuMu CC COH","f");
  tLMode_T_SelCoh->AddEntry(hT_SignalELS_SelCoh,"NuMu CC NonSig","f");
  tLMode_T_SelCoh->AddEntry(hT_OtherNuNC_SelCoh,"Other Nu NC","f");
  tLMode_T_SelCoh->AddEntry(hT_Cosmic_SelCoh,"InEvent Cosmic","f");
  tLMode_T_SelCoh->AddEntry(hT_InTime_SelCoh,"InTime Cosmic","f");

  THStack *hStackMode_T_SelCoh = new THStack("hStackMode_T_SelCoh", "");
  hStackMode_T_SelCoh->Add(hT_SignalCOH_SelCoh);
  hStackMode_T_SelCoh->Add(hT_SignalQEL_SelCoh);
  hStackMode_T_SelCoh->Add(hT_SignalMEC_SelCoh);
  hStackMode_T_SelCoh->Add(hT_SignalRES_SelCoh);
  hStackMode_T_SelCoh->Add(hT_SignalDIS_SelCoh);
  hStackMode_T_SelCoh->Add(hT_SignalELS_SelCoh);
  hStackMode_T_SelCoh->Add(hT_OtherNuNC_SelCoh);
  hStackMode_T_SelCoh->Add(hT_InTime_SelCoh);
  hStackMode_T_SelCoh->Add(hT_Cosmic_SelCoh);

  new TCanvas;
  hT_All_SelCoh->Draw("hist");
  hStackMode_T_SelCoh->Draw("hist same");
  hT_All_SelCoh->Draw("hist same");
  tLMode_T_SelCoh->Draw();
  gPad->Print("Spectrum_RecoT_Modes_CoherentStudy_SelCoh.pdf");

  // zoom
  TLegend *tLMode_Tzoom_SelCoh = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_Tzoom_SelCoh->SetHeader("|t| in COH Selection","C");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalQEL_SelCoh,"NuMu CC QEL","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalMEC_SelCoh,"NuMu CC MEC","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalRES_SelCoh,"NuMu CC RES","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalDIS_SelCoh,"NuMu CC DIS","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalCOH_SelCoh,"NuMu CC COH","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_SignalELS_SelCoh,"NuMu CC NonSig","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_OtherNuNC_SelCoh,"Other Nu NC","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_Cosmic_SelCoh,"InEvent Cosmic","f");
  tLMode_Tzoom_SelCoh->AddEntry(hTzoom_InTime_SelCoh,"InTime Cosmic","f");

  THStack *hStackMode_Tzoom_SelCoh = new THStack("hStackMode_Tzoom_SelCoh", "");
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalCOH_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalQEL_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalMEC_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalRES_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalDIS_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_SignalELS_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_OtherNuNC_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_InTime_SelCoh);
  hStackMode_Tzoom_SelCoh->Add(hTzoom_Cosmic_SelCoh);

  new TCanvas;
  hTzoom_All_SelCoh->Draw("hist");
  hStackMode_Tzoom_SelCoh->Draw("hist same");
  hTzoom_All_SelCoh->Draw("hist same");
  tLMode_Tzoom_SelCoh->Draw();
  gPad->Print("Spectrum_RecoTzoom_Modes_CoherentStudy_SelCoh.pdf");

  // --> with cheated proton cut at 40 MeV
  TLegend *tLMode_Tzoom_ChtCoh = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_Tzoom_ChtCoh->SetHeader("|t| in COH Selection","C");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalQEL_ChtCoh,"NuMu CC QEL","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalMEC_ChtCoh,"NuMu CC MEC","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalRES_ChtCoh,"NuMu CC RES","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalDIS_ChtCoh,"NuMu CC DIS","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalCOH_ChtCoh,"NuMu CC COH","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_SignalELS_ChtCoh,"NuMu CC NonSig","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_OtherNuNC_ChtCoh,"Other Nu NC","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_Cosmic_ChtCoh,"InEvent Cosmic","f");
  tLMode_Tzoom_ChtCoh->AddEntry(hTzoom_InTime_ChtCoh,"InTime Cosmic","f");

  THStack *hStackMode_Tzoom_ChtCoh = new THStack("hStackMode_Tzoom_ChtCoh", "");
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalCOH_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalQEL_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalMEC_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalRES_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalDIS_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_SignalELS_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_OtherNuNC_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_InTime_ChtCoh);
  hStackMode_Tzoom_ChtCoh->Add(hTzoom_Cosmic_ChtCoh);
  
  new TCanvas;
  hTzoom_All_ChtCoh->Draw("hist");
  hStackMode_Tzoom_ChtCoh->Draw("hist same");
  hTzoom_All_ChtCoh->Draw("hist same");
  tLMode_Tzoom_ChtCoh->Draw();
  gPad->Print("Spectrum_RecoTzoom_Modes_CoherentStudy_ChtCoh.pdf");

  // --> with cheated proton cut at 15 MeV
  TLegend *tLMode_Tzoom_CLECoh = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_Tzoom_CLECoh->SetHeader("|t| in COH Selection","C");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalQEL_CLECoh,"NuMu CC QEL","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalMEC_CLECoh,"NuMu CC MEC","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalRES_CLECoh,"NuMu CC RES","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalDIS_CLECoh,"NuMu CC DIS","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalCOH_CLECoh,"NuMu CC COH","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_SignalELS_CLECoh,"NuMu CC NonSig","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_OtherNuNC_CLECoh,"Other Nu NC","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_Cosmic_CLECoh,"InEvent Cosmic","f");
  tLMode_Tzoom_CLECoh->AddEntry(hTzoom_InTime_CLECoh,"InTime Cosmic","f");

  THStack *hStackMode_Tzoom_CLECoh = new THStack("hStackMode_Tzoom_CLECoh", "");
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalCOH_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalQEL_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalMEC_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalRES_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalDIS_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_SignalELS_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_OtherNuNC_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_InTime_CLECoh);
  hStackMode_Tzoom_CLECoh->Add(hTzoom_Cosmic_CLECoh);
  
  new TCanvas;
  hTzoom_All_CLECoh->Draw("hist");
  hStackMode_Tzoom_CLECoh->Draw("hist same");
  hTzoom_All_CLECoh->Draw("hist same");
  tLMode_Tzoom_CLECoh->Draw();
  gPad->Print("Spectrum_RecoTzoom_Modes_CoherentStudy_CLECoh.pdf");

  // RES Selection
  double totalCounts_SelRes = hAll_SelRes->Integral();
  TLegend *tLMode_SelRes = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_SelRes->SetHeader(TString::Format("COH Select = %.1f",hSignalRES_SelRes->Integral()), "C");
  tLMode_SelRes->AddEntry(hSignalQEL_SelRes,TString::Format("NuMu CC QEL: %.1f%%",100.*hSignalQEL_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hSignalMEC_SelRes,TString::Format("NuMu CC MEC: %.1f%%",100.*hSignalMEC_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hSignalRES_SelRes,TString::Format("NuMu CC RES: %.1f%%",100.*hSignalRES_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hSignalDIS_SelRes,TString::Format("NuMu CC DIS: %.1f%%",100.*hSignalDIS_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hSignalCOH_SelRes,TString::Format("NuMu CC COH: %.1f%%",100.*hSignalCOH_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hSignalELS_SelRes,TString::Format("NuMu CC NonSig: %.1f%%",100.*hSignalELS_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hOtherNuNC_SelRes,TString::Format("Other Nu NC: %.1f%%",100.*hOtherNuNC_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hCosmic_SelRes,TString::Format("InEvent Cosmic: %.1f%%",100.*hCosmic_SelRes->Integral()/totalCounts_SelRes),"f");
  tLMode_SelRes->AddEntry(hInTime_SelRes,TString::Format("InTime Cosmic: %.1f%%",100.*hInTime_SelRes->Integral()/totalCounts_SelRes),"f");

  THStack *hStackMode_SelRes = new THStack("hStackMode_SelRes", "");
  hStackMode_SelRes->Add(hSignalCOH_SelRes);
  hStackMode_SelRes->Add(hSignalQEL_SelRes);
  hStackMode_SelRes->Add(hSignalMEC_SelRes);
  hStackMode_SelRes->Add(hSignalRES_SelRes);
  hStackMode_SelRes->Add(hSignalDIS_SelRes);
  hStackMode_SelRes->Add(hSignalELS_SelRes);
  hStackMode_SelRes->Add(hOtherNuNC_SelRes);
  hStackMode_SelRes->Add(hInTime_SelRes);
  hStackMode_SelRes->Add(hCosmic_SelRes);

  new TCanvas;
  hAll_SelRes->Draw("hist");
  hStackMode_SelRes->Draw("hist same");
  hAll_SelRes->Draw("hist same");
  tLMode_SelRes->Draw();
  gPad->Print("Spectrum_Modes_RESStudy_SelRes.pdf");

  // 1MuNP (sometimes called 1Mu1P in the vars...) broken into different classes to understand the background...
  double totalCounts_1MuNP = hSignalwProton_recoMuMom_BaseSliceSel->Integral();
  TLegend *tLMode_1MuNP = new TLegend(0.6,0.6,0.87,0.87);
  tLMode_1MuNP->SetHeader("1MuNP where N!=1 with |p| > 400 MeV/c", "C");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_TrueSignal,"True Signal","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NumuMuWrong,"NumuCC but #mu not signal","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NumuPWrong,"NumuCC but p not signal","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NumuBothWrong,"NumuCC but both wrong","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NumuNonFidVol,"NumuCC but non fiducial","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NumuOther,"NumuCC Other","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_NuOther,"Nu Other","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_Cosmic,"In-spill Cosmic","f");
  tLMode_1MuNP->AddEntry(hSignalwProton_recoMuMom_Cosmic_InTime,"In-time Cosmic","f");

  THStack *hStackMode_1MuNP = new THStack("hStackMode_1MuNP", "");
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_TrueSignal);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NumuMuWrong);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NumuPWrong);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NumuBothWrong);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NumuNonFidVol);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NumuOther);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_NuOther);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_Cosmic);
  hStackMode_1MuNP->Add(hSignalwProton_recoMuMom_Cosmic_InTime);

  new TCanvas;
  hSignalwProton_recoMuMom_BaseSliceSel->Draw("hist");
  hStackMode_1MuNP->Draw("hist same");
  hSignalwProton_recoMuMom_BaseSliceSel->Draw("hist same");
  tLMode_1MuNP->Draw();
  gPad->Print("Spectrum_1Mu1PSelection.pdf");

  std::cout << "" << std::endl;
  std::cout << "Done." << std::endl;
}
