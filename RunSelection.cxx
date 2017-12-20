#include <experimental/optional>
#include <thread>
#include <mutex>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "TSUtil.h"
#include "TSSelection.h"

using namespace std;
using namespace art;

// copied from stack overflow
inline bool isInteger(const std::string & s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char *p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}

bool isNumber(const std::string& s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

    char *end = 0;
    strtod(s.c_str(), &end);
    return end != s.c_str();
}

struct ThreadData {
   int n_selections;
   int n_trials;
   std::vector<int> dataset_id;
   std::vector<float> track_energy_distortion;
   bool track_energy_distortion_by_percent;
   std::vector<float> shower_energy_distortion;
   bool shower_energy_distortion_by_percent;
   bool accept_np;
   bool accept_ntrack;
   std::vector<string> input_files;
   std::vector<TFile *> output_files;     
   int thread_index;
};

void execute_thread(const ThreadData &d, mutex &writer_lock) {
  std::vector<galleryfmwk::TSSelection> selections(d.n_selections);

  writer_lock.lock();
  cout << "Initialize" << endl;
  writer_lock.unlock();
  for (int i = 0; i < d.n_selections; i ++) {
    selections[i].setOutputFile(d.output_files[i]); 
    vector<string> this_input_files(d.input_files);
    selections[i].initialize(this_input_files, &writer_lock);
  }

  writer_lock.lock();
  cout << "Setup" << endl;
  writer_lock.unlock();

  for (int i = 0; i < d.n_selections; i ++) {
    selections[i].setFluxWeightProducer("eventweight");
    selections[i].setEventWeightProducer("mcweight");
    selections[i].setMCTruthProducer("generator");
    selections[i].setMCTrackProducer("mcreco");
    selections[i].setMCShowerProducer("mcreco");
    selections[i].setVerbose(true);
    selections[i].setNTrials(d.n_trials);
    selections[i].setAcceptP(d.accept_np, 2);
    selections[i].setAcceptNTrk(d.accept_ntrack);
   
    if (d.dataset_id.size() > 0) {
        selections[i].setDatasetID( d.dataset_id[i] );
    } 
    if (d.track_energy_distortion.size() > 0) {
      selections[i].setTrackEnergyResolution(d.track_energy_distortion[i], d.track_energy_distortion_by_percent);
    } 
    if (d.shower_energy_distortion.size() > 0) {
      selections[i].setShowerEnergyResolution(d.shower_energy_distortion[i], d.shower_energy_distortion_by_percent);
    }
  }

  writer_lock.lock();
  cout << "Analyze" << endl;
  writer_lock.unlock();
  for (gallery::Event ev(d.input_files) ; !ev.atEnd(); ev.next()) { 
    for (galleryfmwk::TSSelection &selection: selections) {
      selection.analyze(&ev);
    }
  }
  
  writer_lock.lock();
  cout << "Finalize" << endl;
  writer_lock.unlock();
  for (galleryfmwk::TSSelection &selection: selections) {
    selection.finalize();
  }
}

// stand alone executable for TSSelection
int main(int argv, char** argc) {
  if (argv < 3) {
    cout << "Must provide at least three input arguments: fout [fin]" << endl;
    return 1;
  }

  // All of the options
  int n_selections = 1;
  int n_trials = 1;
  int n_threads = 1;
  std::vector<int> dataset_id = {};

  std::vector<float> track_energy_distortion = {};
  bool track_energy_distortion_by_percent = false;

  std::vector<float> shower_energy_distortion = {};  
  bool shower_energy_distortion_by_percent = false;

  bool accept_np = false;
  bool accept_ntrack = false;

  cout << "Processing Runtime arguments" << endl;

  std::vector<TFile *> output_files;
  //We have passed the input file as an argument to the function 
  vector<string> input_files;
  for (int i = 1; i < argv; ++i) {
    if (strcmp(argc[i], "-o") == 0 || strcmp(argc[i], "--output_file") == 0) {
      while (i+1 < argv && argc[i+1][0] != '-') {
        i ++;
        output_files.push_back( new TFile(argc[i], "RECREATE") );
        cout << "Output File " << output_files.size() << ": " << argc[i] << endl;
      }
      continue;
    }
    if (strcmp(argc[i], "-n") == 0 || strcmp(argc[i], "--Nselections") == 0) {
      i ++;
      n_selections = stoi(argc[i]);
      cout << "N Selections: " << n_selections << endl;
      continue;
    }
    if (strcmp(argc[i], "-t") == 0 || strcmp(argc[i], "--Ntrials") == 0) {
      i ++;
      n_trials = stoi(argc[i]);
      cout << "N Trials: " << n_trials << endl;
      continue;
    }
    if (strcmp(argc[i], "-T") == 0 || strcmp(argc[i], "--NThreads") == 0) {
      i ++;
      n_threads = stoi(argc[i]);
      cout << "N Threads: " << n_threads << endl;
      continue;
    }
    if (strcmp(argc[i], "-d") == 0 || strcmp(argc[i], "--datasetId") == 0) {
      while ( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        dataset_id.push_back( stoi(argc[i]) );
      }
      cout << "Dataset Id: " ;
      for (auto d: dataset_id) {
          cout << d << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_energy_distortion") == 0) {
      while( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        track_energy_distortion.push_back( stof(argc[i]) );
      }
      cout << "Track Energy Distortion: ";
      for (auto dist: track_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_edist_by_percent") == 0) {
      track_energy_distortion_by_percent = true;
      cout << "Track Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argc[i], "--S_energy_distortion") == 0) {
      while ( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        shower_energy_distortion.push_back( stof(argc[i]) );
      }
      cout << "Shower Energy Distortion: ";
      for (auto dist: shower_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argc[i], "--S_edist_by_percent") == 0) {
      shower_energy_distortion_by_percent = true;
      cout << "Shower Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argc[i], "--accept_np") == 0) {
      accept_np = true;
      cout << "Accept Np" << endl;
      continue;
    }
    if (strcmp(argc[i], "--accept_ntrk") == 0) {
      accept_ntrack = true;
      cout << "Accept Ntrck" << endl;
      continue;
    }
    input_files.push_back(string(argc[i]));
    cout << "Input file: " << argc[i] << endl;
  }

  assert(track_energy_distortion.size() == n_selections
      || track_energy_distortion.size() == 0);
  assert(shower_energy_distortion.size() == n_selections
      || shower_energy_distortion.size() == 0);
  assert(dataset_id.size() == n_selections
      || dataset_id.size() == 0);
  assert(output_files.size() == n_selections * n_threads);
  n_threads = min((int)input_files.size(), n_threads);

  vector<ThreadData> thread_data_list;
  for (int i = 0; i < n_threads; i ++) {
    int n_input_files = input_files.size() / n_threads;
    if (input_files.size() % n_threads != 0) n_input_files ++;

    int min_input_files = n_input_files * i;
    int max_input_files = min(n_input_files*(i+1), (int)input_files.size());

    thread_data_list.push_back({
      n_selections,
      n_trials,
      vector<int>(dataset_id),
      vector<float>(track_energy_distortion),
      track_energy_distortion_by_percent,
      vector<float>(shower_energy_distortion),
      shower_energy_distortion_by_percent,
      accept_np,
      accept_ntrack,
      vector<string>(&input_files[min_input_files], &input_files[max_input_files]),
      vector<TFile *>(&output_files[i*n_selections], &output_files[(i+1)*n_selections]),
      i,
    });
  }

  if (n_threads == 1) {
    mutex writer_lock;
    execute_thread(thread_data_list[0], writer_lock);
  }
  else {
    vector<thread> threads;
    mutex writer_lock;
    for (auto& d: thread_data_list) {
      threads.push_back(thread(execute_thread, ref(d), ref(writer_lock)));  
    }
    for (auto& t: threads) {
      t.join();
    }
  }
}

