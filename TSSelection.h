/**
 * \file TSSelection.h
 * \brief A truth-based 1l1p event selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef GALLERY_FMWK_TSSELECTION_H
#define GALLERY_FMWK_TSSELECTION_H

#include <map>
#include <string>
#include <random>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "Analysis/ana_base.h"
#include "TSUtil.h"

class TDatabasePDG;
class TFile;
class TH2F;
class TNtuple;
class TTree;

namespace galleryfmwk {

/**
 * \class TSSelection
 * \brief Truth-based selection approximating 1l1p
 */
class TSSelection : galleryfmwk::ana_base {

public:
  // A structure to hold temporary track/shower data during processing
  struct PIDParticle {
    int pdg;
    int pdgtrue;
    TLorentzVector p;
    double evis;
    double eccqe;
    double len;
    bool exiting;

    // Output stream operator to print a PIDParticle
    friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
  };

  TSSelection() : _verbose(false) {}

  bool initialize();

  bool analyze(gallery::Event* ev);

  bool finalize();

  void setVerbose(bool b) { _verbose = b; }

  // Set the producers for data products
  void setFluxWeightProducer(std::string s) { _fw_producer = s; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }

  void setShowerEnergyResolution(float);
  void setTrackEnergyResolution(float);
  float nextShowerEnergyDistortion();
  float nextTrackEnergyDistortion();

  // Set a numeric dataset ID, which is written into the tree as a tag
  void setDatasetID(int id) { _dataset_id = id; }

  // Utility function to test if a list of particles is 1l1p
  static bool is1l1p(std::vector<PIDParticle>& p, int lpdg);

  // Apply track cuts
  static inline bool goodTrack(const sim::MCTrack& t, const simb::MCTruth& truth, float energy_distortion=0.) {
    return (!t.empty() &&
            tsutil::isFromNuVertex(truth, t) &&
            t.Process() == "primary" &&
            t.Start().E() - tsutil::get_pdg_mass(t.PdgCode()) + energy_distortion >= 60);
  }

  // Apply shower cuts
  static inline bool goodShower(const sim::MCShower& s, const simb::MCTruth& truth, float energy_distortion=0.) {
    return (tsutil::isFromNuVertex(truth, s) &&
            s.Process() == "primary" &&
            (s.Start().E() - tsutil::get_pdg_mass(s.PdgCode())) + energy_distortion >= 30);
  }

  // A structure used to hold TTree output
  struct OutputData {
    OutputData() {
      weights = NULL;
    }
    int nupdg;
    double enu;
    double q2;
    double w;
    double q0;
    double q3;
    int int_type;
    int int_mode;
    bool ccnc;
    double eccqe;
    double ep;
    int ppdg;
    double elep;
    double thetalep;
    double philep;
    int lpdg;
    int lpid;
    double llen;
    bool lexit;
    double bnbweight;
    int dataset;
    std::map<std::string, std::vector<double> >* weights;
  };

protected:
  // Data product producers
  std::string _track_producer;
  std::string _fw_producer;
  std::string _ew_producer;
  std::string _mct_producer;
  std::string _mcf_producer;
  std::string _mctrk_producer;
  std::string _mcshw_producer;

  // Counters for efficiency and purity calculations
  size_t true_1e1p;
  size_t good_1e1p;
  size_t miss_1e1p;
  size_t true_1m1p;
  size_t good_1m1p;
  size_t miss_1m1p;

  // Optionally set some energy resolution
  float _shower_energy_resolution;
  std::normal_distribution<float> _shower_energy_distribution;
  float _track_energy_resolution;
  std::normal_distribution<float> _track_energy_distribution;
  // random stuff
  std::mt19937 _gen;



  bool _verbose;  //!< Print verbose output
  int _dataset_id;  //!< An arbitrary numeric ID
  OutputData* _data;  //!< Output data
  TTree* _tree;  //!< Output tree

  TNtuple* _truthtree;
  TNtuple* _mectree;

  TFile* _pdf_file;  //!< File containing dE/dx PDFs
  std::map<int, TH2F*> _trackdedxs;  // Track dE/dx distributions
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_TSSELECTION_H

