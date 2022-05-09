#include "sbnana/CAFAna/Prediction/IPrediction.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "OscLib/IOscCalc.h"

#include <cassert>
#include <iostream>

#include "TDirectory.h"
#include "TObjString.h"

// To implement LoadFrom()
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Prediction/PredictionLinFit.h"
#include "sbnana/CAFAna/Prediction/PredictionNoOsc.h"
#include "sbnana/CAFAna/Prediction/PredictionIncDirt.h"
#include "sbnana/CAFAna/Prediction/PredictionSBNExtrap.h"

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IPrediction> LoadFrom<IPrediction>(TDirectory* dir, const std::string& name)
  {
    TObjString* ptag = (TObjString*)dir->GetDirectory(name.c_str())->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "PredictionNoExtrap") return PredictionNoExtrap::LoadFrom(dir, name);

    if(tag == "PredictionInterp") return PredictionInterp::LoadFrom(dir, name);

    //    if(tag == "PredictionLinFit") return PredictionLinFit::LoadFrom(dir, name);

    if(tag == "PredictionNoOsc") return PredictionNoOsc::LoadFrom(dir, name);

    if(tag == "PredictionIncDirt") return PredictionIncDirt::LoadFrom(dir, name);

    if(tag == "PredictionSBNExtrap") return PredictionSBNExtrap::LoadFrom(dir, name);

    std::cerr << "Unknown Prediction type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictUnoscillated() const
  {
    // Default implementation
    osc::NoOscillations noosc;
    return Predict(&noosc);
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictSyst(osc::IOscCalc* calc,
                                    const SystShifts& syst) const
  {
    assert(syst.IsNominal() && "This Prediction doesn't support PredictSyst(). Did you just mean Predict()?");

    // Default implementation: no treatment of systematics
    return Predict(calc);
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictComponentSyst(osc::IOscCalc* calc,
                                             const SystShifts& syst,
                                             Flavors::Flavors_t flav,
                                             Current::Current_t curr,
                                             Sign::Sign_t sign) const
  {
    assert(syst.IsNominal() && "This Prediction doesn't support PredictSyst(). Did you just mean Predict()?");

    // Default implementation: no treatment of systematics
    return PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  void IPrediction::SaveTo(TDirectory* dir, const std::string& name) const
  {
    assert(0 && "Not implemented");
  }
}
