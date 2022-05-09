#pragma once

#include "sbnana/CAFAna/Prediction/IPrediction.h"

#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"

#include "sbnana/CAFAna/Core/Loaders.h"

namespace ana
{
  /// TODO comment
  class PredictionSBNExtrap: public IPrediction
  {
  public:
    // TODO TODO TODO
    /*
    PredictionSBNExtrap(Loaders& loadersND,
                        Loaders& loadersFD,
                        const HistAxis& axis,
                        const SpillCut& spillcut,
                        const Cut& cut,
                        const SystShifts& shift_mc = kNoShift,
                        const Weight& wei_mc = kUnweighted,
                        const SystShifts& shift_data = kNoShift,
                        const Weight& wei_data = kUnweighted);
    */
    virtual ~PredictionSBNExtrap();

    virtual Spectrum Predict(osc::IOscCalc* calc) const override;

    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;

    OscillatableSpectrum ComponentCC(int from, int to) const override;
    //    Spectrum ComponentNC() const override;

    virtual void SaveTo(TDirectory* dir, const std::string& name) const override;
    static std::unique_ptr<PredictionSBNExtrap> LoadFrom(TDirectory* dir, const std::string& name);

    PredictionSBNExtrap() = delete;

  protected:
    PredictionSBNExtrap(const PredictionNoExtrap& pn,
                        const PredictionNoExtrap& pf,
                        const Spectrum& d)
      : fPredND(pn), fPredFD(pf), fDataND(d)
    {
    }

    PredictionNoExtrap fPredND, fPredFD;
    Spectrum fDataND;
  };

  class SBNExtrapGenerator: public IPredictionGenerator
  {
  public:
    SBNExtrapGenerator(Loaders& loaders_nd,
                       const HistAxis& ax,
                       const SpillCut& spillcut,
                       const Cut& cut,
                       const Weight& wei_mc,
                       const SystShifts& shift_data = kNoShift,
                       const Weight& wei_data = kUnweighted);

    std::unique_ptr<IPrediction> Generate(Loaders& loaders_fd,
                                          const SystShifts& shiftMC = kNoShift) const override;

  protected:
    Loaders& fLoadersND;
    HistAxis fAxis;
    SpillCut fSpillCut;
    Cut fCut;
    Var fWeightMC;
    SystShifts fShiftData;
    Var fWeightData;
  };
}
