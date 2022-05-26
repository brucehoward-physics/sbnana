#include "cafanacore/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  class FluxTimesNuclei: public Spectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    FluxTimesNuclei(INuTruthSource& src, const Binning& bins,
                    const NuTruthCut& fidvol, int pdg, const NuTruthWeight& wgt = kNuTruthUnweighted);

    TH1D* ToTH1(double pot,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);
  protected:
    int fPdg;
  };

  class EnsembleFluxTimesNuclei: public EnsembleSpectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    EnsembleFluxTimesNuclei(INuTruthEnsembleSource& src, const Binning& bins,
                            const NuTruthCut& fidvol, int pdg, const NuTruthWeight& wgt = kNuTruthUnweighted);

    TH1D* ToTH1(double pot,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);

    void SaveTo(TDirectory* dir, const std::string& name) const;

    static std::unique_ptr<EnsembleFluxTimesNuclei> LoadFrom(TDirectory* dir, const std::string& name);

    /// Convert an EnsembleSpectrum (i.e. a \ref EnsembleFluxTimesNuclei) into an ensemble where every
    /// bin within a given universe is the integral of the EnsembleSpectrum
    EnsembleSpectrum MakeTotalFlux(const HistAxis& ax) const;

  protected:
    /// Helper for LoadFrom()
    EnsembleFluxTimesNuclei(const EnsembleSpectrum* spec, const int pdg);

    int fPdg;
  };


}
