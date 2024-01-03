#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/SBNAna/Vars/NuMIXSecTruthVars.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnanaobj/StandardRecord/SRSlice.h"

namespace ana
{
  /// bool to determine if object is in fiducial volume
  bool isInFV (double x, double y, double z);
  bool isInAV (double x, double y, double z);

  /// bool to determine if object is in containment volume
  bool isContainedVol (double x, double y, double z);

  /// Utilities for PFParticle loops
  bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk );
  bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk );
  bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw );
  bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk );

  bool IsValidTrkIdx( const caf::SRSlice* slice, const unsigned int idxTrk );
  bool IsTracklikeTrack( const caf::SRSlice* slice, const unsigned int idxTrk );
  bool IsShowerlike( const caf::SRSlice* slice, const unsigned int idxShw );
  bool IsPrimaryPFP( const caf::SRSlice* slice, const unsigned int idxTrk );

  /// \ref Var that is a dummy var that returns 1 for "IsFHC" or 0 for "!IsFHC" for example
  extern const Var kNuMIDummyVar1;
  extern const Var kNuMIDummyVar0;

  //// Utilities for chi2
  double GetChi2MIP(const caf::Proxy<caf::SRTrackCalo>& calo);

  /// \ref SpillVar for trigger time (check if the implementation is only comaptible for emulated trigger and fix if so...)
  extern const SpillVar kNuMISpillTriggerTime;

  /// \ref Var for the muon candidate index
  extern const Var kNuMIMuonCandidateIdx;

  /// \ref Var for the proton candidate index
  extern const Var kNuMIProtonCandidateIdx;

  /// \ref MultiVar for the charged pion candidate index
  extern const MultiVar kNuMIChargedPionCandidateIdxs;

  /// \ref MultiVar for the proton candidate indices
  extern const MultiVar kNuMIPhotonCandidateIdxs;

  /// kinematic/output variables
  // Muon momentum
  extern const Var kNuMIMuonCandidateRecoP;
  extern const Var kNuMIMuonTrueP;
  // Muon length
  extern const Var kNuMIRecoMuonLength;
  extern const Var kNuMITrueMuonLength;
  // Muon transverse momentum
  extern const Var kNuMIRecoMuonPt;
  extern const Var kNuMITrueMuonPt;
  // Proton momentum
  extern const Var kNuMIProtonCandidateRecoP;
  extern const Var kNuMIProtonTrueP;
  // Proton transverse momentum
  extern const Var kNuMIRecoProtonPt;
  extern const Var kNuMITrueProtonPt;
  // Proton legnth
  extern const Var kNuMIRecoProtonLength;
  extern const Var kNuMITrueProtonLength;
  // Muon angle w.r.t. beam
  extern const Var kNuMIRecoCosThBeam;
  extern const Var kNuMITrueCosThBeam;
  // Muon angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  extern const Var kNuMIRecoCosThVtx;
  extern const Var kNuMITrueCosThVtx;
  // Proton angle w.r.t. beam
  extern const Var kNuMIProtonRecoCosThBeam;
  extern const Var kNuMIProtonTrueCosThBeam;
  // Proton angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  extern const Var kNuMIProtonRecoCosThVtx;
  extern const Var kNuMIProtonTrueCosThVtx;
  // Angle btw muon and proton
  extern const Var kNuMIRecoCosThMuP;
  extern const Var kNuMITrueCosThMuP;

  // TKI variables
  // - delta PT
  extern const Var kNuMIRecodeltaPT;
  extern const Var kNuMITruedeltaPT;
  // - delta PTx
  extern const Var kNuMIRecodeltaPTx;
  extern const Var kNuMITruedeltaPTx;
  // - delta PTy
  extern const Var kNuMIRecodeltaPTy;
  extern const Var kNuMITruedeltaPTy;
  // - delta alphaT
  extern const Var kNuMIRecodeltaalphaT;
  extern const Var kNuMITruedeltaalphaT;
  // - delta phiT
  extern const Var kNuMIRecodeltaphiT;
  extern const Var kNuMITruedeltaphiT;

  // Sideband vars: pi0
  extern const Var kNuMILeadingPhotonCandidateE;
  extern const Var kNuMILeadingPhotonCandidateTrueE;
  extern const Var kNuMISecondaryPhotonCandidateE;
  extern const Var kNuMISecondaryPhotonCandidateTrueE;
  extern const Var kNuMIPhotonCandidatesOpeningAngle;

  // Sideband vars: pi+-
  extern const Var kNuMILeadingChargedPionCandidateInd;
  extern const Var kNuMILeadingChargedPionCandidateLength;
  extern const Var kNuMILeadingChargedPionCandidateNDaughter;
  extern const Var kNuMILeadingChargedPionCandidateMatchedPDG;
  extern const Var kNuMILeadingChargedPionCandidateNCollectionHit;
  extern const Var kNuMILeadingChargedPionCandidateMIPChi2;

  // Investigate shower vars for the cut
  extern const MultiVar kNuMIShowerCandidateIdxs;
  extern const Var kNuMILeadingShowerCandidateLen;
  extern const Var kNuMILeadingShowerCandidateGap;
  extern const Var kNuMILeadingShowerCandidateEColl;
  extern const Var kNuMILeadingShowerCandidateNHitsColl;
  extern const Var kNuMILeadingShowerCandidateOpenAngle;
  extern const Var kNuMILeadingShowerCandidateTrkFitEColl;
}
