#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "cafanacore/Progress.h"
#include "cafanacore/ReweightableSpectrum.h"
#include "cafanacore/SAMProjectSource.h"
#include "cafanacore/Spectrum.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

namespace ana
{
  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::string& wildcard, int max)
    : SpectrumLoaderBase(wildcard), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::vector<std::string>& fnames, int max)
    : SpectrumLoaderBase(fnames), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader()
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader* SpectrumLoader::FromSAMProject(const std::string& proj,
                                                 int fileLimit)
  {
    SpectrumLoader* ret = new SpectrumLoader;
    ret->fWildcard = "project "+proj;
    ret->fFileSource = std::unique_ptr<IFileSource>(new SAMProjectSource(proj, fileLimit));
    return ret;
  }

  //----------------------------------------------------------------------
  SpectrumLoader::~SpectrumLoader()
  {
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::Go()
  {
    if(fGone){
      std::cerr << "Error: can only call Go() once on a SpectrumLoader" << std::endl;
      abort();
    }
    fGone = true;

    const int Nfiles = NFiles();

    Progress* prog = 0;

    caf::SRBranchRegistry::clear();

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      if(Nfiles >= 0 && !prog){
        unsigned int totsinks = 0;
        for(const ISpillSink* s: fSinks) totsinks += s->NSinks();

        prog = new Progress(TString::Format("Filling %u spectra from %d files matching '%s'", totsinks, Nfiles, fWildcard.c_str()).Data());
      }

      HandleFile(f, Nfiles == 1 ? prog : 0);

      if(Nfiles > 1 && prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end for fileIdx

    if(prog){
      prog->Done();
      delete prog;
    }

    StoreExposures(); // also triggers the POT printout

    for(ISpillSink* sink: fSinks) sink->RemoveSource(this);
    fSinks.clear();
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::PrintGraph(std::ostream& os) const
  {
    os << "digraph{" << std::endl;
    os << "comment = \"Render me with a command like: dot -Tpdf graph.dot > graph.pdf\"" << std::endl << std::endl;
    beta::Passthrough<caf::SRSpillProxy>::PrintGraph(os);
    os << "}";
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::PrintGraph() const
  {
    PrintGraph(std::cout);
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());

    // Test for multi (trees are inside a directory) or single (tree is at top
    // level) tree cases.
    TDirectory* dir = 0;
    TTree* tr = 0;

    TObject* obj = f->Get("recTree");
    assert(obj); // Must exist in one form or the other

    // It might seem like you could use GetObject() to do the type-checking
    // directly, but that method seems to have a memory leak in the case that
    // the types don't match.
    if(obj->ClassName() == std::string("TTree")){
      tr = (TTree*)obj;
    }
    else{
      dir = (TDirectory*)obj;
      tr = (TTree*)dir->Get("rec");
      assert(tr);
    }

    // We try to access this field for every record. It was only added to the
    // files in late 2021, and we don't want to render all earlier files
    // unusable at a stroke. This logic can safely be removed once all extant
    // files have such a field (estimate mid-2022?)
    const bool has_husk = tr->GetLeaf("rec.hdr.husk");

    const caf::CAFType type = caf::GetCAFType(dir, tr);

    long n;
    caf::SRSpillProxy sr(dir, tr, "rec", n, 0);

    //    FloatingExceptionOnNaN fpnan;

    long Nentries = tr->GetEntries();
    if (max_entries != 0 && max_entries < Nentries) Nentries = max_entries;

    for(n = 0; n < Nentries; ++n){
      if(type != caf::kFlatMultiTree) tr->LoadTree(n); // for all single-tree modes

      if(sr.hdr.first_in_subrun){
        fPOT += sr.hdr.pot;
        // TODO think about if this should be gated behind first_in_file. At
        // the moment I think these will be synonymous. And despite the comment
        // on hdr.pot, I think it may be file-based in practice too.
        fNGenEvt += sr.hdr.ngenevt;
      }

      // This record was only kept as a receptacle for exposure information. It
      // shouldn't be included in any selected spectra.
      if(has_husk && sr.hdr.husk) continue;

      HandleRecord(&sr, 1);

      if(prog) prog->SetProgress(double(n)/Nentries);
    } // end for n
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    if(fabs(fPOT - fPOTFromHist)/std::min(fPOT, fPOTFromHist) > 0.001){
      std::cout << fPOT << " POT from hdr differs from " << fPOTFromHist << " POT from the TotalPOT histogram!" << std::endl;
      abort();
    }

    std::cout << fPOT << " POT over " << fNGenEvt << " readouts" << std::endl;

    FillPOT(fPOT);
    FillLivetime(fNGenEvt);
  }

} // namespace
