#include <TClonesArray.h>      // or TClonesArray
#include <TGenericClassInfo.h> // for TGenericClassInfo
#include <TMath.h>             // for Sqrt
#include <TRandom.h>           // for TRandom, gRandom
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <algorithm> // std::sort
#include <vector>    // std::vector
#include <cmath>
#include "FairMCEventHeader.h"    // for FairMCEventHeader
#include "FairLink.h"             // for FairLink
#include "FairRunAna.h"           // for FairRunAna
#include "FairRootManager.h"      // for FairRootManager
#include "ConvDriftTubeRawData.h" // for Conversion
#include "DriftTube.h"            // for Drift Tube detector

// FIXME cleanup headers

ConvDriftTubeRawData::ConvDriftTubeRawData()
   : FairTask("ConvDriftTubeRawData"), fSNDTree(nullptr), fMiniDTTree(nullptr), fDigiDriftTube(nullptr),
     MiniDTeventNumber(0), MatchedEntries(0)
{
}

ConvDriftTubeRawData::~ConvDriftTubeRawData() {}

InitStatus ConvDriftTubeRawData::Init()
{
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      LOG(FATAL) << "ConvDriftTubeRawData::Init: RootManager not instantiated!";
   }

   // Get the DriftTube detector from the list of globals
   DriftTubeDet = dynamic_cast<DriftTube *>(gROOT->GetListOfGlobals()->FindObject("DriftTube"));

   // Input raw data file is read from the FairRootManager
   // This allows to have it in custom format, e.g. have arbitary names of TTrees

   // // raw SND
   // TFile *f0 = dynamic_cast<TFile *>(ioman->GetObject("data"));
   // fSNDTree = (TTree *)f0->Get("data");
   // fSNDTree->GetBranch("evt_timestamp");

   // converted SND
   TFile *f0 = dynamic_cast<TFile *>(ioman->GetObject("rawConv"));
   fSNDTree = (TTree *)f0->Get("rawConv");
   // fSNDTree->GetBranch("EventHeader.fEventTime");

   // auto fMiniDT = static_cast<TFile *>(TFile::Open("/afs/cern.ch/user/g/guiducci/public/snd/mdt_tree_norb100_novl50_run_010743.root"));
   auto fMiniDT = static_cast<TFile *>(TFile::Open("/eos/user/g/guiducci/temp-analysis/results_run_010988/mdt_tree_Norb10_Novl50_run_010988.root"));
   // fMiniDTTree = static_cast<TTree *>(fMiniDT->Get("minidt_hits"));
   fMiniDTTree = static_cast<TTree *>(fMiniDT->Get("minidt"));

   // Register the output
   fDigiDriftTube = new TClonesArray("DriftTubeHit");
   ioman->Register("Digi_DriftTubeHits", "DigiDriftTubeHit_det", fDigiDriftTube, kTRUE);

   // Get the FairLogger
   FairLogger *logger = FairLogger::GetLogger();

   eventNumber = fnStart;

   return kSUCCESS;
}

void ConvDriftTubeRawData::Exec(Option_t * /*opt*/)
{

   fDigiDriftTube->Clear("C"); //

   // Delete pointer map elements
   for (auto it : digiDTStore) {
      delete it.second;
   }
   digiDTStore.clear();

   // run the conversion
   Process();
   eventNumber++;
}

void ConvDriftTubeRawData::Process()
{
   int indexDriftTube{}; // index of DT hits
   int detID;
   fSNDTree->GetEvent(eventNumber);
   auto eventTimestamp = fSNDTree->GetLeaf("EventHeader.fEventTime")->GetValue();

   TTreeReader MiniDTReader(fMiniDTTree);

   // TTreeReaderValue<int> n_hits(MiniDTReader, "n_hits");
   TTreeReaderArray<long long> hit_orbit(MiniDTReader, "hit_orbit");
   // TTreeReaderArray<int> hit_bx(MiniDTReader, "hit_bx");
   // TTreeReaderArray<int> hit_tdc(MiniDTReader, "hit_tdc");
   TTreeReaderArray<double> hit_timestamp(MiniDTReader, "hit_timestamp");
   // TTreeReaderArray<int> hit_chamber(MiniDTReader, "hit_chamber");
   // TTreeReaderArray<int> hit_layer(MiniDTReader, "hit_layer");
   // TTreeReaderArray<int> hit_wire(MiniDTReader, "hit_wire");
   // TTreeReaderValue<int> n_tpgs(MiniDTReader, "n_tpgs");
   // TTreeReaderArray<int> tpg_t0(MiniDTReader, "tpg_t0");
   // TTreeReaderArray<int> tpg_position(MiniDTReader, "tpg_position");
   // TTreeReaderArray<int> tpg_slope(MiniDTReader, "tpg_slope");
   // TTreeReaderArray<int> tpg_chi2(MiniDTReader, "tpg_chi2");
   // TTreeReaderArray<int> tpg_hit_l1_time(MiniDTReader, "tpg_hit_l1_time");
   // TTreeReaderArray<int> tpg_hit_l1_lat(MiniDTReader, "tpg_hit_l1_lat");
   // TTreeReaderArray<int> tpg_hit_l1_wire(MiniDTReader, "tpg_hit_l1_wire"); // layer is inferred form the hit order in
   // the tpg index[0, 4] --> layer[0, 4] TTreeReaderArray<int> tpg_hit_l1_valid(MiniDTReader, "tpg_hit_l1_valid");
   // TTreeReaderArray<int> tpg_hit_l2_time(MiniDTReader, "tpg_hit_l2_time");
   // TTreeReaderArray<int> tpg_hit_l2_lat(MiniDTReader, "tpg_hit_l2_lat");
   // TTreeReaderArray<int> tpg_hit_l2_wire(MiniDTReader, "tpg_hit_l2_wire");
   // TTreeReaderArray<int> tpg_hit_l2_valid(MiniDTReader, "tpg_hit_l2_valid");
   // TTreeReaderArray<int> tpg_hit_l3_time(MiniDTReader, "tpg_hit_l3_time");
   // TTreeReaderArray<int> tpg_hit_l3_lat(MiniDTReader, "tpg_hit_l3_lat");
   // TTreeReaderArray<int> tpg_hit_l3_wire(MiniDTReader, "tpg_hit_l3_wire");
   // TTreeReaderArray<int> tpg_hit_l3_valid(MiniDTReader, "tpg_hit_l3_valid");
   // TTreeReaderArray<int> tpg_hit_l4_time(MiniDTReader, "tpg_hit_l4_time");
   // TTreeReaderArray<int> tpg_hit_l4_lat(MiniDTReader, "tpg_hit_l4_lat");
   // TTreeReaderArray<int> tpg_hit_l4_wire(MiniDTReader, "tpg_hit_l4_wire");
   // TTreeReaderArray<int> tpg_hit_l4_valid(MiniDTReader, "tpg_hit_l4_valid");
   // TTreeReaderArray<int> tpg_chamber(MiniDTReader, "tpg_chamber");
   // TTreeReaderArray<int> tpg_orbit(MiniDTReader, "tpg_orbit");
   TTreeReaderArray<double> tpg_timestamp(MiniDTReader, "tpg_timestamp");

   MiniDTReader.SetEntry(MiniDTeventNumber);
   while (MiniDTReader.Next()) {
      if ((std::count_if(hit_orbit.begin(), hit_orbit.end(), [](int n) { return n < 0; }) == 0) && (hit_timestamp.GetSize() > 0)) {
         TTreeReaderArray<double>::iterator HitTimestampIt = std::min_element(hit_timestamp.begin(), hit_timestamp.end());
         // double MinHitTimestamp = *HitTimestampIt;
         double t_diff_hit = *HitTimestampIt - static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));
         // double t_diff_tpg = tpg_timestamp[0] - static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));

         // if ((abs(t_diff_hit) < 4e-6) || abs(t_diff_tpg < 1.6e-6)) { // così abbiamo finestra prima e dopo hit/tpg
         // MiniDT
         if (abs(t_diff_hit) < 6e-6) { // così abbiamo finestra prima e dopo hit/tpg MiniDT
            std::cout << "t_diff_hit: " << t_diff_hit << '\n';
            // std::cout << "t_diff_tpg: " << t_diff_tpg << '\n';
            // std::cout << "t SND: " << eventTimestamp / (4 * 40.0789 * 1e6) << " t MiniDT hit: " << hit_timestamp[0]
            // << " t MiniDT tpg: " << tpg_timestamp[0] << '\n';
            std::cout << "t SND: " << eventTimestamp / (4 * 40.0789 * 1e6) << " t MiniDT hit: " << *HitTimestampIt << '\n';
            digiDTStore[detID] = new DriftTubeHit(detID);
            MiniDTeventNumber = MiniDTReader.GetCurrentEntry(); // tenere questo indice come nuovo start per reader
            ++MatchedEntries;
            break;
         } else {
            continue;
         }
      } else {
         break;
      }
   }

   // // Loop over hits per event!
   // for (int n = 0; n < fSNDTree->GetLeaf("nHits")->GetValue(); n++) { // FIXME - this is simply an example

   //    // FIXME add conversion below

   //    // store the hits internally per event
   //    if (digiDTStore.count(detID) == 0) {
   //       // make the hit via the hit constructor
   //       // e.g. digiDTStore[detID] = new DriftTubeHit(detID,xxx);
   //    }
   //    // would be nice to have a link btw raw hit and converted one - below is example
   //    // from the SciFi where we use the DAQ RO params
   //    // digiDTStore[detID]->SetDaqID(sipm_number,n, board_id, tofpet_id, tofpet_channel);
   // } // end loop over hits in the event

   // // write the hits to the output
   // for (auto it_detID : digiDTStore) {
   //    (*fDigiDriftTube)[indexDriftTube] = digiDTStore[it_detID.first];
   //    indexDriftTube += 1;
   // }
   // LOG(INFO) << fnStart + 1 << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
   
   // LOG(INFO) << eventNumber << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
   UpdateInput(eventNumber);
}

void ConvDriftTubeRawData::UpdateInput(int NewStart)
{
   fSNDTree->Refresh();
   eventNumber = NewStart;
}

void ConvDriftTubeRawData::PrintMatchedEntries()
{
   std::cout << MatchedEntries << " matched entries out of " << fMiniDTTree->GetEntries()
             << " MiniDT Tree total entries. \n";
}