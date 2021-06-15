// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SSHiggs : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SSHiggs);

	// lambda function to sort by pT
	ParticleSorter ptsorter =  [](const Particle& a, const Particle& b) {
		return a.pt() > b.pt();
	};


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      // lepton pT
      book(_h["l0_pt"], "l0_pt", 200, 0.0, 1000.0);
      book(_h["l1_pt"], "l1_pt", 200, 0.0, 1000.0);

      // MET
      book(_h["MET"], "MET", 100, 0.0, 500.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons(ptsorter);

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV);

      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      // Veto event if there are no b-jets
      //if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      //if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;
      const auto MET = apply<MissingMomentum>(event, "MET").missingPt();

      // === fill histograms ===
      _h["l0_pt"]->fill(leptons[0].pT()/GeV);
      _h["l1_pt"]->fill(leptons[1].pT()/GeV);

      _h["MET"]->fill(MET/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
	
		for (auto hist : _h) {
			scale(hist.second, crossSection() / sumOfWeights());
		}

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    //map<string, Profile1DPtr> _p;
    //map<string, CounterPtr> _c;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(SSHiggs);

}
