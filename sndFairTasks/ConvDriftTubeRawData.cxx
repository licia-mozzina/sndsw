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

   // converted SND
   TFile *f0 = dynamic_cast<TFile *>(ioman->GetObject("rawConv"));
   fSNDTree = (TTree *)f0->Get("rawConv");
   // fSNDTree->GetBranch("EventHeader.fEventTime");

   fMiniDTChain = (TChain *)ioman->GetObject("MiniDTChain");

   // fMiniDTChain = new TChain("minidt_hits"); // FIXME delete?
   // fMiniDTChain->Add("/eos/user/g/guiducci/temp-analysis/after_ts1_analysis/minidt_run_011833_trees/minidt_run_011833_hits.root");
   // std::array<int, 68> nFiles;
   // std::iota(nFiles.begin(), nFiles.end(), 1);
   // for (auto i : nFiles) {
   //    if (i > 0) {
   //       char fMiniDT[200];
   //       sprintf(fMiniDT, "/eos/user/g/guiducci/temp-analysis/after_ts1_analysis/minidt_run_011833_trees/minidt_run_011833_hits_%i.root", i);
   //       fMiniDTChain->Add(fMiniDT);
   //    } else continue;
   // }

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

   TTreeReaderValue<double> hit_timestamp(MiniDTReader, "hit_timestamp");
   TTreeReaderValue<int> hit_chamber(MiniDTReader, "hit_chamber");
   TTreeReaderValue<int> hit_layer(MiniDTReader, "hit_layer");
   TTreeReaderValue<int> hit_wire(MiniDTReader, "hit_wire");

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
         // std::cout << MatchedHits << '\n';
      } else if ((*hit_timestamp - SNDtimestamp) > 1e-6) {
         for (auto it_detID : digiDTStore) {
            (*fDigiDriftTube)[indexDriftTube] = digiDTStore[it_detID.first];
            indexDriftTube += 1;
         }
         if (digiDTStore.size() != 0) {
            // std::cout << digiDTStore.size() << " matched hits!\n";
         }
         // MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
         MatchedHits = 0;
         ++MatchedEntries;
         break;
      } else {
         continue;
      }
   }

   // LOG(INFO) << eventNumber << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
   UpdateInput(eventNumber);
}

void ConvDriftTubeRawData::UpdateInput(int NewStart)
{
//   fSNDTree->Refresh();
   eventNumber = NewStart;
}

