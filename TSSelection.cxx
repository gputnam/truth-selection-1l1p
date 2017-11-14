#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// includes for random draws from gaussian
#include <random>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "gallery/Event.h"

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "TSUtil.h"
#include "TSSelection.h"

namespace galleryfmwk {
void TSSelection::setShowerEnergyResolution(float res) {
  _shower_energy_resolution = res;
  _shower_energy_distribution = std::normal_distribution<float>(0., res);
}

void TSSelection::setTrackEnergyResolution(float res) {
  _track_energy_resolution = res;
  _track_energy_distribution = std::normal_distribution<float>(0., res);
}

float TSSelection::nextTrackEnergyDistortion() {
  if (_track_energy_resolution < 1e-4)
    return 0.;
  return _track_energy_distribution( _gen );
}

float TSSelection::nextShowerEnergyDistortion() {
  if (_shower_energy_resolution < 1e-4)
    return 0.;
  return _shower_energy_distribution( _gen );
}

bool TSSelection::is1l1p(std::vector<PIDParticle>& p, int lpdg) {
  // Are there exactly two particles ID'ed?
  if (p.size() > 2) {
    return false;
  }

  // Count protons and the chosen lepton type
  size_t np = 0;
  size_t nl = 0;
  for (size_t i=0; i<p.size(); i++) {
    if (p[i].pdg == 2212) {
      np++;
    }
    if (p[i].pdg == lpdg) {
      nl++;
    }
  }

  return np == 1 && nl == 1;
}


bool TSSelection::initialize() {
  // Load track dE/dx distributions from file
  _pdf_file = TFile::Open("./data/dedx_pdfs.root");
  assert(_pdf_file->IsOpen());

  TIter next(_pdf_file->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    const char* name = key->GetName();
    TObjArray* tokens = TString(name).Tokenize("_");
    TString htype = ((TObjString*)tokens->At(0))->GetString();
    int pdg = atoi(((TObjString*)tokens->At(1))->GetString());

    if (htype.Contains("htrackdedx")) {
      TH2F* h = (TH2F*) _pdf_file->Get(name);
      if (h->Integral() == 0 || pdg < 0 || pdg > 10000) { continue; }

      // Ignore a few low bins
      for (int i=0; i<h->GetNbinsX(); i++) {
        for (int j=0; j<h->GetNbinsY(); j++) {
          if (i < 2 || j < 2) {
            h->SetBinContent(i, j, 0);
          }
        }
      }
      _trackdedxs[pdg] = h;
    }
  }

  // Initialize dataset identifier
  _dataset_id = -1;

  // Initialize event counters
  true_1e1p = 0;
  good_1e1p = 0;
  miss_1e1p = 0;
  true_1m1p = 0;
  good_1m1p = 0;
  miss_1m1p = 0;

  // initialize resolutions to 0 (perfect resolution)
  _shower_energy_resolution = 0.;
  _track_energy_resolution = 0.;
  _shower_energy_distribution = std::normal_distribution<float>(0.0, 0.0);
  _track_energy_distribution = std::normal_distribution<float>(0.0, 0.0);

  // setting up random # stuff
  std::random_device rd;
  _gen = std::mt19937( rd() );;

  // Set up the output trees
  assert(_fout);
  _fout->cd();

  _data = new OutputData;
  _tree = new TTree("data", "");
  _tree->Branch("nupdg", &_data->nupdg);
  _tree->Branch("enu", &_data->enu);
  _tree->Branch("q2", &_data->q2);
  _tree->Branch("w", &_data->w);
  _tree->Branch("q0", &_data->q0);
  _tree->Branch("q3", &_data->q3);
  _tree->Branch("int", &_data->int_type);
  _tree->Branch("mode", &_data->int_mode);
  _tree->Branch("ccnc", &_data->ccnc);
  _tree->Branch("eccqe", &_data->eccqe);
  _tree->Branch("ep", &_data->ep);
  _tree->Branch("ppdg", &_data->ppdg);
  _tree->Branch("elep", &_data->elep);
  _tree->Branch("thetalep", &_data->thetalep);
  _tree->Branch("philep", &_data->philep);
  _tree->Branch("lpdg", &_data->lpdg);
  _tree->Branch("lpid", &_data->lpid);
  _tree->Branch("llen", &_data->llen);
  _tree->Branch("lexit", &_data->lexit);
  _tree->Branch("bnbweight", &_data->bnbweight);
  _tree->Branch("dataset", &_data->dataset);
  _tree->Branch("weights", &_data->weights);

  _truthtree = new TNtuple("truth", "", "nupdg:enu:ccnc:int:mode:w:q2:lpdg:elep:tlep:npip:npim:npi0:np:nn:fw:ttrk:rtrk:texit:tshr:rshr:sexit");
  _mectree = new TNtuple("mec", "", "nupdg:enu:ccnc:mode:w:q2:lpdg:tlep:ep0:ep1:ep2:ep3:ep4");

  return true;
}


bool TSSelection::analyze(gallery::Event* ev) {
  // Get handles for event data
  art::InputTag gtruth_tag(_mct_producer);
  auto const& gtruth_list = \
    (*ev->getValidHandle<std::vector<simb::GTruth> >(gtruth_tag));

  art::InputTag eventweight_tag(_ew_producer);
  auto const& eventweights_list = \
    (*ev->getValidHandle<std::vector<evwgh::MCEventWeight> >(eventweight_tag));

  art::InputTag mctruth_tag(_mct_producer);
  auto const& mctruth_list = \
    (*ev->getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(_mcshw_producer);
  auto const& mcshower_list = \
    (*ev->getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(_mctrk_producer);
  auto const& mctrack_list = \
    (*ev->getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  // Sanity checks
  assert(_fout);
  assert(mctruth_list.size() == gtruth_list.size());

  // BNB flux weight
  double wbnb = 1.0;
  if (!eventweights_list.empty() &&
      eventweights_list[0].fWeight.find("bnbcorrection_FluxHist") != eventweights_list[0].fWeight.end()) {
    wbnb = eventweights_list[0].fWeight.at("bnbcorrection_FluxHist")[0];
  }

  _fout->cd();

  // Loop through MC truth interactions
  for(size_t i=0; i<mctruth_list.size(); i++) {
    const simb::MCTruth& mctruth = mctruth_list.at(i);
    const simb::GTruth& gtruth = gtruth_list.at(i);

    size_t ntracks = 0, nshowers = 0;

    // Keep track of event particle content (currently a little redundant)
    std::vector<PIDParticle> particles_found;
    std::vector<PIDParticle> particles_true;

    // Get vertex-associated contained tracks
    for (size_t j=0; j<mctrack_list.size(); j++) {
      const sim::MCTrack& mct = mctrack_list.at(j);

      float energy_distortion = nextTrackEnergyDistortion();

      // Apply track cuts
      if (!goodTrack(mct, mctruth, energy_distortion)) {
        continue;
      }

      // Track length
      double s = 0;
      TLorentzVector pos = mct.End().Position();
      for (long k=mct.size()-2; k>=0; k--) {
        s += (pos.Vect() - mct[k].Position().Vect()).Mag();
        pos = mct[k].Position();
      }

      // don't apply energy distortion to the "true" particle data
      particles_true.push_back({
        mct.PdgCode(),
        mct.PdgCode(),
        mct.Start().Momentum(),
        mct.Start().E() - tsutil::get_pdg_mass(mct.PdgCode()),
        tsutil::eccqe(mct.Start().Momentum()),
        s,
        !tsutil::inFV(mct)
      });

      ntracks++;

      // KS test on dE/dx vs. range to pick the best match for PID
      TH2F* htemp = tsutil::HistDEdx(mct);
      double ks_best = 0;
      int pdg_best = -999;

      for (auto const& it : _trackdedxs) {
        if (htemp->Integral() == 0 || it.second->Integral() == 0) { continue; }

        double ks = it.second->KolmogorovTest(htemp);
        if (ks > ks_best) {
          ks_best = ks;
          pdg_best = it.first;
        }
      }

      delete htemp;

      // Un-PID "protons" that are too long or short
      if (pdg_best == 2212 && (s > 80 || s < 12)) {
        pdg_best = -888;
      }

      // Call all remaining unmatched tracks protons
      if (pdg_best == -999) {
        pdg_best = 2212;
      }

      particles_found.push_back({
        pdg_best,
        mct.PdgCode(),
        mct.Start().Momentum(),
        mct.Start().E() - tsutil::get_pdg_mass(mct.PdgCode()) + energy_distortion,
        tsutil::eccqe(mct.Start().Momentum()),
        s,
        !tsutil::inFV(mct)
      });
    }


    // Get vertex-associated contained showers
    for (size_t j=0; j<mcshower_list.size(); j++) {
      const sim::MCShower& mcs = mcshower_list.at(j);

      float energy_distortion = nextShowerEnergyDistortion();

      // Apply shower cuts
      if (!goodShower(mcs, mctruth, energy_distortion)) {
        continue;
      }

      // don't apply energy distortion to the "true" particle data
      particles_true.push_back({
        mcs.PdgCode(),
        mcs.PdgCode(),
        mcs.Start().Momentum(),
        mcs.Start().E() - tsutil::get_pdg_mass(mcs.PdgCode()),
        tsutil::eccqe(mcs.Start().Momentum()),
        -1,
        !tsutil::inFV(mcs)
      });

      nshowers++;

      // Guess the PDG based on shower dE/dx (cut at 3.5)
      int pdg_best = (mcs.dEdx() < 3.5 ? 11 : 22);

      particles_found.push_back({
        pdg_best,
        mcs.PdgCode(),
        mcs.Start().Momentum(),
        // @ANDY: IS THIS A BUG???
        // previous lines have this energy as:
        // mcs.Start().E() - tsutil::get_pdg_mass(mcs.PdgCode())
        mcs.Start().E() + energy_distortion,
        tsutil::eccqe(mcs.Start().Momentum()),
        -1,
        !tsutil::inFV(mcs)
      });
    }

    // Classify the event (found/true 1l1p/1m1p)
    // "True" 1l1p here means there are one true l and one true p that pass the
    // track/shower cuts (i.e. are in within this specific signal definition).
    bool f_1e1p = is1l1p(particles_found, 11);
    bool t_1e1p = is1l1p(particles_true, 11);
    bool f_1m1p = is1l1p(particles_found, 13);
    bool t_1m1p = is1l1p(particles_true, 13);

    // Where have all the muons gone?
    //if (t_1m1p && !f_1m1p) {
    //  std::cout << "1m1p missed" << std::endl;
    //  std::cout << "true: n=" << particles_true.size() << ": ";
    //  for (size_t i=0; i<particles_true.size(); i++) {
    //    std::cout << particles_true[i].pdg << "/" << particles_true[i].pdgtrue << "(" << particles_true[i].evis << ") ";
    //  }
    //  std::cout << std::endl;

    //  std::cout << "found: n=" << particles_found.size() << ": ";
    //  for (size_t i=0; i<particles_found.size(); i++) {
    //    std::cout << particles_found[i].pdg << "/" << particles_found[i].pdgtrue << "(" << particles_found[i].evis << ") ";
    //  }
    //  std::cout << std::endl;
    //}

    if (t_1e1p) true_1e1p++;
    if (f_1e1p && t_1e1p) good_1e1p++;
    if (f_1e1p && !t_1e1p) miss_1e1p++;

    if (t_1m1p) true_1m1p++;
    if (f_1m1p && t_1m1p) good_1m1p++;
    if (f_1m1p && !t_1m1p) miss_1m1p++;

    // Print out PID information mis-IDs
    if ((f_1e1p && !t_1e1p) || (f_1m1p && !t_1m1p)) {
      std::cout << "true: " << mctruth.GetNeutrino().Nu().E() * 1000
                << "[" << mctruth.GetNeutrino().InteractionType() << "] ";
      for (size_t k=0; k<particles_true.size(); k++) {
        std::cout << particles_true[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "est: " << ntracks << " tracks, "
                           << nshowers << " showers; ";
      for (size_t k=0; k<particles_found.size(); k++) {
        std::cout << particles_found[k] << " ";
      }
      std::cout << std::endl;
    }

    // Write event to output tree for found 1l1p events
    if (f_1e1p || f_1m1p) {
      double eccqe=-1, ep=-1, ppdg=-1, elep=-1, thetalep=-1, philep=-1, lpdg=-1, lpid=-1, llen=-1, lexit=-1;

      for (size_t k=0; k<particles_found.size(); k++) {
        if (particles_found[k].pdg == 2212) {
          ep = particles_found[k].evis;
          ppdg = particles_found[k].pdgtrue;
        }
        else {
          eccqe = particles_found[k].eccqe;
          elep = particles_found[k].evis;
          thetalep = particles_found[k].p.Theta();
          philep = particles_found[k].p.Phi();
          lpid = particles_found[k].pdg;
          lpdg = particles_found[k].pdgtrue;
          llen = particles_found[k].len;
          lexit = particles_found[k].exiting;
        }
      }

      const simb::MCNeutrino& nu = mctruth.GetNeutrino();
      const simb::MCParticle& pnu = nu.Nu();
      const simb::MCParticle& plep = nu.Lepton();
      TLorentzVector xp = (pnu.Momentum() - plep.Momentum());

      std::map<std::string, std::vector<double> > wgh;
      if (!eventweights_list.empty()) {
        wgh = eventweights_list[0].fWeight;
      }

      _data->nupdg = nu.Nu().PdgCode();
      _data->enu = nu.Nu().E();
      _data->q2 = nu.QSqr();
      _data->w = nu.W();
      _data->q0 = xp.E();
      _data->q3 = xp.Vect().Mag();
      _data->int_type = nu.InteractionType();
      _data->int_mode = nu.Mode();
      _data->ccnc = nu.CCNC();
      _data->eccqe = eccqe;
      _data->ep = ep;
      _data->ppdg = ppdg;
      _data->elep = elep;
      _data->thetalep = thetalep;
      _data->philep = philep;
      _data->lpdg = lpdg;
      _data->lpid = lpid;
      _data->llen = llen;
      _data->lexit = lexit;
      _data->bnbweight = wbnb;
      _data->dataset = _dataset_id;
      _data->weights = &wgh;
      _tree->Fill();
    }

    // Fill the event truth tree
    float vtt[22] = {
      (float) mctruth.GetNeutrino().Nu().PdgCode(),
      (float) mctruth.GetNeutrino().Nu().E(),
      (float) mctruth.GetNeutrino().CCNC(),
      (float) mctruth.GetNeutrino().InteractionType(),
      (float) mctruth.GetNeutrino().Mode(),
      (float) mctruth.GetNeutrino().W(),
      (float) mctruth.GetNeutrino().QSqr(),
      (float) mctruth.GetNeutrino().Lepton().PdgCode(),
      (float) mctruth.GetNeutrino().Lepton().E(),
      (float) gtruth.fGint,
      (float) gtruth.fNumPiPlus,
      (float) gtruth.fNumPiMinus,
      (float) gtruth.fNumPi0,
      (float) gtruth.fNumProton,
      (float) gtruth.fNumNeutron,
      (float) wbnb,
      (float) mctrack_list.size(),
      (float) ntracks,
      (float) 0,
      (float) mcshower_list.size(),
      (float) nshowers,
      (float) 0
    };
    _truthtree->Fill(vtt);

    // Fill tree for MEC events with sorted proton energies
    if (mctruth.GetNeutrino().Mode() == simb::kMEC) {
      std::vector<float> epmec(5, -1);
      for (size_t k=0; k<mctrack_list.size(); k++) {
        const sim::MCTrack& t = mctrack_list[k];
        if (t.PdgCode() == 2212 && t.Process() == "primary") {
          double ke = t.Start().E() - tsutil::get_pdg_mass(t.PdgCode());
          epmec.push_back(ke);
        }
      }

      std::sort(epmec.begin(), epmec.end(), std::greater<>());

      float v2[13] = {
        (float) mctruth.GetNeutrino().Nu().PdgCode(),
        (float) mctruth.GetNeutrino().Nu().E(),
        (float) mctruth.GetNeutrino().CCNC(),
        (float) mctruth.GetNeutrino().Mode(),
        (float) mctruth.GetNeutrino().W(),
        (float) mctruth.GetNeutrino().QSqr(),
        (float) mctruth.GetNeutrino().Lepton().PdgCode(),
        (float) mctruth.GetNeutrino().Lepton().E(),
        (float) epmec[0],
        (float) epmec[1],
        (float) epmec[2],
        (float) epmec[3],
        (float) epmec[4]
      };
      _mectree->Fill(v2);
    }
  }

  return true;
}


bool TSSelection::finalize() {
  // Print out statistics
  std::cout << "1e1p true: " << true_1e1p
            << ", good: " << good_1e1p
            << ", miss: " << miss_1e1p
            << std::endl;

  std::cout << "1e1p eff: "
            << 1.0 * (good_1e1p + miss_1e1p) / true_1e1p
            << std::endl;

  std::cout << "1e1p pur: "
            << 1.0 * good_1e1p / (good_1e1p + miss_1e1p)
            << std::endl;

  std::cout << "1m1p true: " << true_1m1p
            << ", good: " << good_1m1p
            << ", miss: " << miss_1m1p << std::endl;

  std::cout << "1m1p eff: "
            << 1.0 * (good_1m1p + miss_1m1p) / true_1m1p
            << std::endl;

  std::cout << "1m1p pur: "
            << 1.0 * good_1m1p / (good_1m1p + miss_1m1p)
            << std::endl;

  // Write output to ROOT file
  if (_fout) {
    _fout->cd();
    _tree->Write();
    _truthtree->Write();
    _mectree->Write();
  }

  return true;
}


std::ostream& operator<<(std::ostream& os, const TSSelection::PIDParticle& dt) {
  os << dt.pdg << "(" << dt.evis << ")";
  return os;
}

}  // namespace gallery_fmwk

