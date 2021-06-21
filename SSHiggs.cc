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

      IdentifiedFinalState electrons(fs);
      electrons.acceptIdPair(PID::ELECTRON);
      declare(electrons, "Electrons");
      
      IdentifiedFinalState muons(fs);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "Muons");

      DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), true, false);
      declare(dressed_muons, "DressedMuons");

      DressedLeptons dressed_electrons(photons, electrons, 0.1, Cuts::open(), true, false);
      declare(dressed_electrons, "DressedElectrons");


      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      // lepton pT
      book(_h["l0_pt"], "l0_pt", 200, 0.0, 1000.0);
      book(_h["l1_pt"], "l1_pt", 200, 0.0, 1000.0);

      book(_h["el0_pt"], "el0_pt", 200, 0.0, 1000.0);
      book(_h["el1_pt"], "el1_pt", 200, 0.0, 1000.0);
      
      book(_h["mu0_pt"], "mu0_pt", 200, 0.0, 1000.0);
      book(_h["mu1_pt"], "mu1_pt", 200, 0.0, 1000.0);

      // MET
      book(_h["MET"], "MET", 100, 0.0, 500.0);

      // object multiplicity
      book(_h["Njet"], "n_j", 11, -0.5, 10.5);
      book(_h["Nlepton"], "n_l", 11, -0.5, 10.5);
      book(_h["Nmuon"], "n_mu", 11, -0.5, 10.5);
      book(_h["Nelectron"], "n_e", 11, -0.5, 10.5);

      // mass plots

      // mass of lead two positive charge leptons
      book(_h["m_pos_leptons"], "m_pos_leptons", 200, 0., 1000.);
      // mass of lead two negative charge leptons
      book(_h["m_neg_leptons"], "m_neg_leptons", 200, 0., 1000.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // === retrieve objects ===
      
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

      std::vector<DressedLepton> muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons(ptsorter);
      std::vector<DressedLepton> electrons = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons(ptsorter);

      // === event selection ===

      // Veto event if there are no b-jets
      //if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      //if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;
      const auto MET = apply<MissingMomentum>(event, "MET").missingPt();

      // object multiplicities
      const int nLeptons   = leptons.size();
      const int nMuons     = muons.size();
      const int nElectrons = electrons.size();
      const int nJets      = jets.size();

      // calculate di-lepton masses (for four-lepton events)
      double m_pos_leptons = -999;
      double m_neg_leptons = -999;

      // first, sort leptons into positive and negative vectors
      std::vector<DressedLepton> positiveLeptons;
      std::vector<DressedLepton> negativeLeptons;
      for (auto lep : leptons) {
         if (lep.charge() > 0) {
            positiveLeptons.push_back(lep);
         }
         if (lep.charge() < 0) {
            negativeLeptons.push_back(lep);
         }
      }

      // calculate dilepton mass, if we have at least two leptons
      if (positiveLeptons.size() >= 2) {
         m_pos_leptons = ( positiveLeptons[0].momentum() + positiveLeptons[1].momentum() ).mass();
      }

      if (negativeLeptons.size() >= 2) {
         m_neg_leptons = ( negativeLeptons[0].momentum() + negativeLeptons[1].momentum() ).mass();
      }

      // === fill histograms ===
      if (nLeptons > 0) {
         _h["l0_pt"]->fill(leptons[0].pT()/GeV);
      }

      if (nLeptons > 1) {
         _h["l1_pt"]->fill(leptons[1].pT()/GeV);
      }

      if (nElectrons > 0) {
         _h["el0_pt"]->fill(electrons[0].pT()/GeV);
      }

      if (nElectrons > 1) {
         _h["el1_pt"]->fill(electrons[1].pT()/GeV);
      }

      if (nMuons > 0) {
         _h["mu0_pt"]->fill(muons[0].pT()/GeV);
      }

      if (nMuons > 1) {
         _h["mu1_pt"]->fill(muons[1].pT()/GeV);
      }

      _h["MET"]->fill(MET/GeV);

      _h["Njet"]->fill(nJets);
      _h["Nlepton"]->fill(nLeptons);
      _h["Nmuon"]->fill(nMuons);
      _h["Nelectron"]->fill(nElectrons);

      _h["m_pos_leptons"]->fill(m_pos_leptons);
      _h["m_neg_leptons"]->fill(m_neg_leptons);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
	
		for (auto hist : _h) {
			scale(hist.second, crossSection() / sumOfWeights());
		}

		// write to std::cout so that we can capture xsec info
		std::cout << "Rivet cross section: " << crossSection() << std::endl;

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
