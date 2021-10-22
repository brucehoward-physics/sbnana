#pragma once

#include "sbnana/CAFAna/StandardRecord/Proxy/FwdDeclare.h"

namespace ana{

  /// Whether this is a primary particle or generated by a secondary interaction
  bool IsPrimary(const caf::SRTrueParticleProxy& p);

  /// Whether this particle should have a bragg peak in the detector
  bool HasBraggPeak(const caf::SRTrueParticleProxy& p);

  /// Whether this particle was generated by genie (as opposed to geant or corsika)
  bool IsGenie(const caf::SRTrueParticleProxy& p);

  /// Whether this is a stable particle as generated by genie
  bool IsStable(const caf::SRTrueParticleProxy& p);

}
