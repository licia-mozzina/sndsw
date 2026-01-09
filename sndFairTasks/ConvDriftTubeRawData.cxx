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
#include <array>     // std::array
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

#define HCELL 1.3
#define WCELL 4.2

ConvDriftTubeRawData::ConvDriftTubeRawData()
   : FairTask("ConvDriftTubeRawData"), fSNDTree(nullptr), fMiniDTChain(nullptr), fDigiDriftTube(nullptr),
     MiniDTeventNumber(0)
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

   // // Delete pointer map elements
   // for (auto it : digiDTStore) {
   //    delete it.second;
   // }
   // digiDTStore.clear();

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

   TTreeReader MiniDTReader(fMiniDTChain);

   TTreeReaderValue<double> hit_timestamp(MiniDTReader, "hit_timestamp");
   TTreeReaderValue<int> hit_chamber(MiniDTReader, "hit_chamber");
   TTreeReaderValue<int> hit_layer(MiniDTReader, "hit_layer");
   TTreeReaderValue<int> hit_wire(MiniDTReader, "hit_wire");

   MiniDTReader.SetEntry(MiniDTeventNumber);
   int MatchedHits {};
   while (MiniDTReader.Next()) {
      auto SNDtimestamp = static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));
      if (((*hit_timestamp - SNDtimestamp) > -2e-7) && ((*hit_timestamp - SNDtimestamp) < 1e-6)) {
         if (MatchedHits == 0) {
            MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
         }
         // if (digiDTStore.count(MatchedHits) == 0) {                                             // so we do not kill subsequent events in same cell
         //    detID = SetDetID(*hit_chamber, *hit_layer, *hit_wire);
         //    digiDTStore[MatchedHits] = new DriftTubeHit(detID, *hit_timestamp - SNDtimestamp);
         // }
         detID = SetDetID(*hit_chamber, *hit_layer, *hit_wire);
         (*fDigiDriftTube)[MatchedHits] = new DriftTubeHit(detID, *hit_timestamp - SNDtimestamp);
         ++MatchedHits;
      } else if ((*hit_timestamp - SNDtimestamp) > 1e-6) {
         // for (auto const& it_detID : digiDTStore) { // si può eliminare?
         //    (*fDigiDriftTube)[indexDriftTube] = it_detID.second; // tanto non abbiamo un matching 1-1 con hit e detID, visto che possono esserci più hit nella stessa cella ed è solo riempito in maniera sequenziale e non ci interessa matching con indice di fdigidrifttube
         //    indexDriftTube += 1;
         // }
         // if (digiDTStore.size() != 0) {
         //    std::cout << digiDTStore.size() << " matched hits!\n";
         // }
         
         // run laterality computation and hit redefinition 
        std::vector<std::vector<int>> Clusters = FindClusters(fDigiDriftTube);   // separate clusters for each chamber!
         // for (clusters) {
            // FindLateralityHough(cluster, fDigiDriftTube);
         // }
         MatchedHits = 0;
         break;
      } else {
         continue;
      }
   }

   LOG(INFO) << eventNumber << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
   UpdateInput(eventNumber);
}

// void ConvDriftTubeRawData::Process()
// {
//    int indexDriftTube{}; // index of DT hits
//    int detID;            // assegna valore!!
//    fSNDTree->GetEvent(eventNumber);
//    auto eventTimestamp = fSNDTree->GetLeaf("EventHeader.fEventTime")->GetValue();

//    TTreeReader MiniDTReader(fMiniDTChain);

//    TTreeReaderValue<double> hit_timestamp(MiniDTReader, "hit_timestamp");
//    TTreeReaderValue<int> hit_chamber(MiniDTReader, "hit_chamber");
//    TTreeReaderValue<int> hit_layer(MiniDTReader, "hit_layer");
//    TTreeReaderValue<int> hit_wire(MiniDTReader, "hit_wire");

//    MiniDTReader.SetEntry(MiniDTeventNumber);
//    int MatchedHits = 0;
//    while (MiniDTReader.Next()) {
//       auto SNDtimestamp = static_cast<double>(eventTimestamp / (4 * 40.0789 * 1e6));
//       if (((*hit_timestamp - SNDtimestamp) > -2e-7) && ((*hit_timestamp - SNDtimestamp) < 1e-6)) {
//          if (MatchedHits == 0) {
//             MiniDTeventNumber = MiniDTReader.GetCurrentEntry();
//          }
//          if (digiDTStore.count(MatchedHits) == 0) {
//             digiDTStore[MatchedHits] = new DriftTubeHit(MatchedHits, *hit_timestamp - SNDtimestamp, *hit_chamber, *hit_layer, *hit_wire);
//          }
//          ++MatchedHits;
//          std::cout << MatchedHits << '\n';
//       } else if ((*hit_timestamp - SNDtimestamp) > 1e-6) {
//          for (auto it_detID : digiDTStore) {
//             (*fDigiDriftTube)[indexDriftTube] = digiDTStore[it_detID.first];
//             indexDriftTube += 1;
//          }
//          if (digiDTStore.size() != 0) {
//             std::cout << digiDTStore.size() << " matched hits!\n";
//          }
//          MatchedHits = 0;
//          break;
//       } else {
//          continue;
//       }
//    }

//    LOG(INFO) << eventNumber << " events processed out of " << fSNDTree->GetEntries() << " number of events in file.";
//    UpdateInput(eventNumber);
// }


void ConvDriftTubeRawData::UpdateInput(int NewStart)
{
//   fSNDTree->Refresh();
   eventNumber = NewStart;
}

int ConvDriftTubeRawData::SetDetID(const int& chamber, const int& layer, const int& wire) 
{
   int detID = 40000 + 1000 * (1 - chamber) + 100 * (3 - layer) + (15 - wire);
   return detID;
}

struct ConvDriftTubeRawData::HitPoint {
   int L; 
   int C;
};

std::vector<ConvDriftTubeRawData::HitPoint> ConvDriftTubeRawData::GetNeighbours(int L, int C) {
   std::vector<HitPoint> neighbours;
   neighbours.reserve(10);

   neighbours.push_back({L, C - 1}); 
   neighbours.push_back({L, C + 1});

   for (int nextL : {L -1, L + 1}) {
      if (nextL < 0 || nextL >= 4) continue;

      neighbours.push_back({nextL, C});

      if (L % 2 == 0) {
         neighbours.push_back({nextL, C + 1});
      } else {
         neighbours.push_back({nextL, C - 1});
      }
   }

   for (int nextL : {L - 2, L + 2}) {
      if (nextL < 0 || nextL >= 4) continue;

      neighbours.push_back({nextL, C});
      neighbours.push_back({nextL, C - 1});
      neighbours.push_back({nextL, C + 1});
   }

   return neighbours;
}

std::vector<std::vector<int>> ConvDriftTubeRawData::FindClusters(const TClonesArray * hits) {
   std::array<int, 2 * 4 * 16> grid;
   grid.fill(-1);

   auto gridIdx = [&](int p, int l, int c) -> int& {return grid[p * (4 * 16) + l * 16 + c];};

   int nHits {static_cast<int>(hits->GetEntries())};
   if (nHits == 0) return {};
   std::vector<bool> visited(nHits, false);
   std::vector<std::vector<int>> allClusters;

   // Fill grid 
   for (int i = 0; i != nHits; ++i) {
      auto hit = dynamic_cast<DriftTubeHit*>(hits->At(i));
      int p = hit->GetPlane();
      int l = hit->GetLayer();
      int c = hit->GetCell();

      if (p >= 0 && p < 2 && l >= 0 && l < 4 && c >= 0 && c < 16) {
         gridIdx(p, l, c) = i;
      }
   }

   // Find cluster
   for (int p = 0; p != 2; ++p) {
      for (int l = 0; l != 4; ++l) {
         for (int c = 0; c != 16; ++c) {
            int startIdx {gridIdx(p, l, c)};

            if (startIdx == -1 || visited[startIdx]) continue;

            std::vector<int> blob;
            std::vector<int> stack = {startIdx};
            visited[startIdx] = true;

            while (!stack.empty()) {
               int currIdx {stack.back()};
               stack.pop_back();
               blob.push_back(currIdx);

               auto hit = dynamic_cast<DriftTubeHit*>(hits->At(currIdx));
               auto neighbours = GetNeighbours(hit->GetLayer(), hit->GetCell());

               for (const auto& nb : neighbours) {
                  if (nb.C >= 0 && nb.C < 16) {
                     int nbIdx {gridIdx(p, nb.L, nb.C)};

                     if (nbIdx != -1 && !visited[nbIdx]) {
                        visited[nbIdx] = true;
                        stack.push_back(nbIdx);
                     }
                  }
               }
            }

            // Split clusters
            if (blob.size() >= 3) {
               std::array<std::vector<int>, 4> layerHits;
               for (int idx : blob) {
                  auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
                  layerHits[hit->GetLayer()].push_back(idx);
               }

               for (auto& layerVec : layerHits) {
                  if (layerVec.empty()) layerVec.push_back(-1);
               }

               for (int i0 : layerHits[0]) {
                  for (int i1 : layerHits[1]) {
                     for (int i2 : layerHits[2]) {
                        for (int i3 : layerHits[3]) {
                           std::vector<int> candidate;

                           candidate.reserve(4);
                           if (i0 != -1) candidate.push_back(i0);
                           if (i1 != -1) candidate.push_back(i1);
                           if (i2 != -1) candidate.push_back(i2);
                           if (i3 != -1) candidate.push_back(i3);

                           if (candidate.size() >= 3) {
                              allClusters.push_back(candidate);
                           }

                        }
                     }
                  }
               }

            }
         }
      }
   }

   for (auto& vector : allClusters) {
      for (auto& v : vector) {
         auto hit = dynamic_cast<DriftTubeHit*>(hits->At(v));
         std::cout << hit->GetPlane() << " " << hit->GetLayer() << " " << hit->GetCell() << "\n";
      }
      std::cout << "Next: \n";
   }
   std::cout << '\n';

   return allClusters;
}

void ConvDriftTubeRawData::FindLateralityHough(const TClonesArray * hits) { // return type may become void

   int nHits {static_cast<int>(hits->GetEntries())};

   if (nHits < 3) return;

   const double SlopeRange {2.};
   const int SlopeBins {100};
   const double dslope {SlopeRange / SlopeBins};

   const double dintercept {0.1};
   const double InterceptOffset = (1.5 * HCELL) * round(SlopeRange / 2);
   const double InterceptBound {(2 * InterceptOffset) + 16.5 * WCELL}; // CHECK AND CHANGE 
   const int InterceptBins {round(InterceptBound / dintercept)};

   
   std::vector<int> Accumulator(SlopeBins * InterceptBins, 0);
   auto bin = [&](int s, int i) -> int& { return Accumulator[s * InterceptBins + i]; };

   for (int i = 0; i != nHits; ++i) {
      auto hit = dynamic_cast<DriftTubeHit*>(hits->At(i)); 
      double y {((hit->GetLayer() - 2) + 0.5) * HCELL};
      double x_offset {(hit->GetCell() + ((hit->GetLayer() % 2 == 1) ? 0.5 : 1.0)) * WCELL};
      for (int lat : {-1, 1}) {
         double x {x_offset + lat * hit->GetTimestamp() * VDRIFT};
         for (int slope_bin = 0; slope_bin != SlopeBins; ++slope_bin) {
            double slope {slope_bin * dslope - round(SlopeRange / 2)};
            double intercept {x + slope * y + InterceptOffset};
            int intercept_bin = round(intercept);
            if (intercept_bin >= 0 && intercept_bin < InterceptBins) {
               bin(slope_bin, intercept_bin)++; 
            }
         }
      } 
   }

   auto it = std::max_element(Accumulator.begin(), Accumulator.end());

   if (*it > 2) {
      size_t max_idx = std::distance(Accumulator.begin(), it);
      auto s_idx = [&](size_t idx) { return idx / InterceptBins; }; // Row
      auto i_idx = [&](size_t idx) { return idx % InterceptBins; }; // Column

      double AccumulatorMaxSlope {s_idx(max_idx) * dslope - round(SlopeRange / 2)};
      double AccumulatorMaxIntercept {i_idx(max_idx) * dintercept - InterceptOffset}; 

      for (int i = 0; i != nHits; ++i) {
         auto hit = dynamic_cast<DriftTubeHit*>(hits->At(i));
         double y {((hit->GetLayer() - 2) + 0.5) * HCELL};
         double x_offset {(hit->GetCell() + ((hit->GetLayer() % 2 == 1) ? 0.5 : 1.0)) * WCELL};
         double x_expected {AccumulatorMaxIntercept - AccumulatorMaxSlope * y};
         double x_drift {hit->GetTimestamp() * VDRIFT};
         if (std::abs(x_expected - x_offset) < WCELL / 2) { // potrebbero esserci cluster da 5 hit
            double residual_left {std::abs(x_expected - (x_offset - x_drift))}; // CHECK WITH SND GEO FOR LEFT RIGHT
            double residual_right {std::abs(x_expected - (x_offset + x_drift))};
            int laterality {(residual_left < residual_right) ? -1 : 1};
            hit->setLaterality(laterality);
         }
      }
   }

   return;
}