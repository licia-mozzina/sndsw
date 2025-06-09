#include <TClonesArray.h>      // or TClonesArray
#include <TGenericClassInfo.h> // for TGenericClassInfo
#include <TMath.h>             // for Sqrt
#include <TRandom.h>           // for TRandom, gRandom
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <algorithm> // std::sort
#include <vector>    // std::vector
#include <cmath>
#include <numeric>   // std::iota
#include <stdio.h>   // sprintf
#include "FairMCEventHeader.h"    // for FairMCEventHeader
#include "FairLink.h"             // for FairLink
#include "FairRunAna.h"           // for FairRunAna
#include "FairRootManager.h"      // for FairRootManager
#include "ConvDriftTubeRawData.h" // for Conversion
#include "DriftTube.h"            // for Drift Tube detector

// FIXME cleanup headers

ConvDriftTubeRawData::ConvDriftTubeRawData()
   : FairTask("ConvDriftTubeRawData"), fSNDTree(nullptr), fMiniDTChain(nullptr), fDigiDriftTube(nullptr),
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

   fMiniDTChain = new TChain("minidt_hits"); // FIXME delete?
   fMiniDTChain->Add("/eos/user/g/guiducci/temp-analysis/trees_minidt_run_011049/minidt_run_011049_hits.root");
   std::array<int, 16> nFiles;
   std::iota(nFiles.begin(), nFiles.end(), 1);
   for (auto i : nFiles) {
      if (i > 0) {
         char fMiniDT[200];
         sprintf(fMiniDT, "/eos/user/g/guiducci/temp-analysis/trees_minidt_run_011049/minidt_run_011049_hits_%i.root", i);
         fMiniDTChain->Add(fMiniDT);
         // fMiniDTChain->Add("/eos/user/g/guiducci/temp-analysis/trees_minidt_run_011049/minidt_run_011049_hits_%i.root");
      } else continue;
   }


   // auto fMiniDT = static_cast<TFile *>(
   //    TFile::Open("/eos/user/g/guiducci/temp-analysis/trees_minidt_run_011049/minidt_run_011049_tpgs_y.root"));
   // fMiniDTTree = static_cast<TTree *>(fMiniDT->Get("minidt_tpgs_y"));

   // Register the output
   fDigiDriftTube = new TClonesArray("DriftTubeHit"); // FIXME dictionary?
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
   int detID;            // assegna valore!!
   fSNDTree->GetEvent(eventNumber);
   auto eventTimestamp = fSNDTree->GetLeaf("EventHeader.fEventTime")->GetValue();

   TTreeReader MiniDTReader(fMiniDTChain);

   // TTreeReaderValue<long long> hit_orbit(MiniDTReader, "hit_orbit");
   // TTreeReaderValue<int> hit_bx(MiniDTReader, "hit_bx");
   // TTreeReaderValue<int> hit_tdc(MiniDTReader, "hit_tdc");
   TTreeReaderValue<double> hit_timestamp(MiniDTReader, "hit_timestamp");
   TTreeReaderValue<int> hit_chamber(MiniDTReader, "hit_chamber");
   TTreeReaderValue<int> hit_layer(MiniDTReader, "hit_layer");
   TTreeReaderValue<int> hit_wire(MiniDTReader, "hit_wire");
   // TTreeReaderValue<int> tpg_t0(MiniDTReader, "tpg_t0");
   // TTreeReaderValue<int> tpg_position(MiniDTReader, "tpg_position");
   // TTreeReaderValue<int> tpg_slope(MiniDTReader, "tpg_slope");
   // TTreeReaderValue<int> tpg_chi2(MiniDTReader, "tpg_chi2");
   // TTreeReaderValue<int> tpg_hit_l1_time(MiniDTReader, "tpg_hit_l1_time");
   // TTreeReaderValue<int> tpg_hit_l1_lat(MiniDTReader, "tpg_hit_l1_lat");
   // TTreeReaderValue<int> tpg_hit_l1_wire(MiniDTReader, "tpg_hit_l1_wire"); // layer is inferred form the hit order in
   // the tpg index[0, 4] --> layer[0, 4] TTreeReaderValue<int> tpg_hit_l1_valid(MiniDTReader, "tpg_hit_l1_valid");
   // TTreeReaderValue<int> tpg_hit_l2_time(MiniDTReader, "tpg_hit_l2_time");
   // TTreeReaderValue<int> tpg_hit_l2_lat(MiniDTReader, "tpg_hit_l2_lat");
   // TTreeReaderValue<int> tpg_hit_l2_wire(MiniDTReader, "tpg_hit_l2_wire");
   // TTreeReaderValue<int> tpg_hit_l2_valid(MiniDTReader, "tpg_hit_l2_valid");
   // TTreeReaderValue<int> tpg_hit_l3_time(MiniDTReader, "tpg_hit_l3_time");
   // TTreeReaderValue<int> tpg_hit_l3_lat(MiniDTReader, "tpg_hit_l3_lat");
   // TTreeReaderValue<int> tpg_hit_l3_wire(MiniDTReader, "tpg_hit_l3_wire");
   // TTreeReaderValue<int> tpg_hit_l3_valid(MiniDTReader, "tpg_hit_l3_valid");
   // TTreeReaderValue<int> tpg_hit_l4_time(MiniDTReader, "tpg_hit_l4_time");
   // TTreeReaderValue<int> tpg_hit_l4_lat(MiniDTReader, "tpg_hit_l4_lat");
   // TTreeReaderValue<int> tpg_hit_l4_wire(MiniDTReader, "tpg_hit_l4_wire");
   // TTreeReaderValue<int> tpg_hit_l4_valid(MiniDTReader, "tpg_hit_l4_valid");
   // TTreeReaderValue<int> tpg_chamber(MiniDTReader, "tpg_chamber");
   // TTreeReaderValue<int> tpg_orbit(MiniDTReader, "tpg_orbit");
   // TTreeReaderValue<double> tpg_timestamp(MiniDTReader, "tpg_timestamp");

   MiniDTReader.SetEntry(MiniDTeventNumber);
   int MatchedHits = 0;
   while (MiniDTReader.Next()) {
      auto SNDtimestamp = static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));
      if (((*hit_timestamp - SNDtimestamp) > -2e-7) && ((*hit_timestamp - SNDtimestamp) < 1e-6)) {
         if (MatchedHits == 0) {
            MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
         }
         if (digiDTStore.count(MatchedHits) == 0) {
            digiDTStore[MatchedHits] = new DriftTubeHit(MatchedHits, *hit_timestamp - SNDtimestamp, *hit_chamber, *hit_layer, *hit_wire);
         }
         ++MatchedHits;
         std::cout << MatchedHits << '\n';
      } else if ((*hit_timestamp - SNDtimestamp) > 1e-6) {
         for (auto it_detID : digiDTStore) {
            (*fDigiDriftTube)[indexDriftTube] = digiDTStore[it_detID.first];
            indexDriftTube += 1;
         }
         if (digiDTStore.size() != 0) {
            std::cout << digiDTStore.size() << " matched hits!\n";
         }
         // MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
         MatchedHits = 0;
         ++MatchedEntries;
         break;
      } else {
         continue;
      }
   }

   // MiniDTReader.SetEntry(MiniDTeventNumber);
   // int MatchedTPGs = 0;
   // while (MiniDTReader.Next()) {
   //    auto SNDtimestamp = static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));
   //    if (((*tpg_timestamp - SNDtimestamp) > -1e-7) && ((*tpg_timestamp - SNDtimestamp) < 3e-7)) {
   //       if (MatchedTPGs == 0) {
   //          MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
   //       }
   //       if (digiDTStore.count(MatchedTPGs) == 0) {
   //          digiDTStore[MatchedTPGs] = new DriftTubeHit(MatchedTPGs, *tpg_timestamp - SNDtimestamp, *tpg_t0, *tpg_position, *tpg_slope, *tpg_chi2);
   //       }
   //       ++MatchedTPGs;
   //       std::cout << MatchedTPGs << '\n';
   //    } else if ((*tpg_timestamp - SNDtimestamp) > 3e-7) {
   //       for (auto it_detID : digiDTStore) {
   //          (*fDigiDriftTube)[indexDriftTube] = digiDTStore[it_detID.first];
   //          indexDriftTube += 1;
   //       }
   //       if (digiDTStore.size() != 0) {
   //          std::cout << digiDTStore.size() << " matched hits!\n";
   //       }
   //       MatchedTPGs = 0;
   //       ++MatchedEntries;
   //       break;
   //    } else {
   //       continue;
   //    }
   // }

   LOG(INFO) << eventNumber << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
   UpdateInput(eventNumber);
}

void ConvDriftTubeRawData::UpdateInput(int NewStart)
{
   fSNDTree->Refresh();
   eventNumber = NewStart;
}

void ConvDriftTubeRawData::PrintMatchedEntries()
{
   std::cout << MatchedEntries << " matched entries out of " << fMiniDTChain->GetEntries()
             << " MiniDT Tree total entries. \n";
}