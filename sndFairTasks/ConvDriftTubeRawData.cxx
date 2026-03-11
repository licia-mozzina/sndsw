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
         detID = SetDetID(*hit_chamber, *hit_layer, *hit_wire);
         (*fDigiDriftTube)[MatchedHits] = new DriftTubeHit(detID, *hit_timestamp - SNDtimestamp);
         ++MatchedHits;
      } else if ((*hit_timestamp - SNDtimestamp) > 1e-6) {  
         // run laterality computation and hit redefinition 
         std::vector<std::vector<int>> hitsClusters = FindClusters(fDigiDriftTube);  
         FindLateralitySlope(fDigiDriftTube, hitsClusters);
         // std::vector<int> latSlope;
         // for (int i = 0; i != fDigiDriftTube->GetEntries(); ++i) {
         //    auto hit = static_cast<DriftTubeHit*>(fDigiDriftTube->At(i));
         //    latSlope.push_back(hit->GetLaterality());
         // }
         // FindLateralityHough(fDigiDriftTube, hitsClusters);
         // for (int i = 0; i != fDigiDriftTube->GetEntries(); ++i) {
         //    auto hit = static_cast<DriftTubeHit*>(fDigiDriftTube->At(i));
         //    std::cout << "latSlope:\t" << latSlope[i] << "latHough:\t" << hit->GetLaterality();
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

   for (int nextL : {L - 1, L + 1}) {
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
                  if (nb.C >= 0 && nb.C < 16) { // bound C, not bounded in GetNeighbours
                     int nbIdx {gridIdx(p, nb.L, nb.C)}; // are we sure we are not mixing hits from diff chambers? yes, the C is bounded and p appears in the function

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

                           if (candidate.size() < 3) continue;

                           bool validChain {true};
                           for (size_t k = 0; k < candidate.size() - 1; ++k) {
                              auto hA = dynamic_cast<DriftTubeHit*>(hits->At(candidate[k]));
                              auto hB = dynamic_cast<DriftTubeHit*>(hits->At(candidate[k + 1]));

                              int dL {hB->GetLayer() - hA->GetLayer()};
                              int dC {hB->GetCell() - hA->GetCell()};

                              if (dL == 1) {
                                 if (hA->GetLayer() % 2 == 0) {
                                    if (dC != 0 && dC != 1) {
                                       validChain = false;
                                       break;
                                    }
                                 } else {
                                    if (dC != 0 && dC != -1) {
                                       validChain = false;
                                       break;
                                    }
                                 }
                              } else if (dL == 2) {
                                 if (std::abs(dC) > 1) {
                                    validChain = false;
                                    break;
                                 }
                              }
                           }

                           if (validChain) {
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

   // for (auto& vector : allClusters) {
   //    for (auto& v : vector) {
   //       auto hit = dynamic_cast<DriftTubeHit*>(hits->At(v));
   //       std::cout << hit->GetPlane() << " " << hit->GetLayer() << " " << hit->GetCell() << "\n";
   //    }
   //    std::cout << "Next: \n";
   // }
   // std::cout << '\n';

   return allClusters;
}


void ConvDriftTubeRawData::FindLateralityHough(const TClonesArray * hits, const std::vector<std::vector<int>>& clusters) { 

   struct TrackResult {
      double slope;
      double intercept;
      int max;
      std::vector<int> indices;
   };

   const auto HCELL = DriftTubeDet->GetConfParF("DriftTube/cellHeight") / 10; // MANCA 1 cm profilo Al
   const auto WCELL = DriftTubeDet->GetConfParF("DriftTube/cellWidth") / 10; 

   const double SlopeRange {2.0};
   const int SlopeBins {100};
   const double dslope {SlopeRange / SlopeBins};

   const double dintercept {0.1};
   const double InterceptOffset = (1.5 * HCELL) * round(SlopeRange / 2);
   const double InterceptBound {(2 * InterceptOffset) + 16.5 * WCELL}; 
   const int InterceptBins {static_cast<int>(InterceptBound / dintercept)};

   static std::vector<int> Accumulator(SlopeBins * InterceptBins, 0);

   std::vector<int> usedBins;
   usedBins.reserve(500);
   std::vector<TrackResult> candidates;
   
   for (auto& hitsIdx : clusters) {
      int nHits {hitsIdx.size()};
      if (nHits < 3) continue;

      TrackResult bestInCluster {0., 0., 0, hitsIdx};

      usedBins.clear();

      for (int idx : hitsIdx) {
         auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx)); 
         double y {((hit->GetLayer() - 2) + 0.5) * HCELL};
         double x_wire {(hit->GetCell() + ((hit->GetLayer() % 2 == 1) ? 0.5 : 1.0)) * WCELL};
         double drift_space {(hit->GetTimestamp() - TPED) * VDRIFT};

         for (int lat : {-1, 1}) {
            double x {x_wire + lat * drift_space};
            if (hit->GetTimestamp() < 0) {
               x = x_wire;
            } else if (hit->GetTimestamp() > WCELL / VDRIFT) {
               x = x_wire + lat * 0.5 * WCELL;
            }

            
            for (int slope_bin = 0; slope_bin != SlopeBins; ++slope_bin) {
               double slope {slope_bin * dslope - round(SlopeRange / 2)};
               double intercept {x + slope * y + InterceptOffset};
               int intercept_bin = round(intercept);
               if (intercept_bin >= 0 && intercept_bin < InterceptBins) {
                  int index {slope_bin * InterceptBins + intercept_bin};
                  if (Accumulator[index] == 0) usedBins.push_back(index);
                  Accumulator[index]++;
               }
            }
         } 
      }

      int max {};
      int bestIdx {-1};

      for (int idx : usedBins) {
         if (Accumulator[idx] > max) {
            max = Accumulator[idx];
            bestIdx = idx;
         }
      }

      // if ((max < 5) && (max > bestInCluster.max)) { // CHECK
      if (max > 2) { // CHECK
         bestInCluster.max = max;
         double maxSlope {(bestIdx / InterceptBins) * dslope - round(SlopeRange / 2)};
         double maxIntercept {(bestIdx / InterceptBins) * dintercept - round(SlopeRange / 2)};
         bestInCluster.slope = maxSlope;
         bestInCluster.intercept = maxIntercept;
         bestInCluster.indices = hitsIdx;
         // std::cout << "results:  " << max << " " << maxSlope << " " << maxIntercept << '\n';
      }
   
      candidates.push_back(bestInCluster);
         
      for (int idx : usedBins) Accumulator[idx] = 0;
      
   }

   std::sort(candidates.begin(), candidates.end(), [](const TrackResult& a, const TrackResult& b) { return a.max > b.max;});

   std::vector<bool> hitUsed(static_cast<int>(hits->GetEntries()), false);

   for (auto& cand : candidates) {
      bool conflict {false};
      
      for (int idx : cand.indices) {
         if (hitUsed[idx]) {
            conflict = true;
            break;
         }
      }

      if (!conflict) {

         for (int idx : cand.indices) {
            auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
            double y {((hit->GetLayer() - 2) + 0.5) * HCELL};
            double x_wire {(hit->GetCell() + ((hit->GetLayer() % 2 == 1) ? 0.5 : 1.0)) * WCELL};
            double drift_space {hit->GetTimestamp() * VDRIFT};

            double x_expected {cand.intercept - cand.slope * y};

            double dist_left {std::abs(x_expected - (x_wire - drift_space))};
            double dist_right {std::abs(x_expected - (x_wire + drift_space))};

            if (dist_left > WCELL / 2) continue;

            int laterality {(dist_left < dist_right) ? -1 : 1}; 

            hit->setLaterality(laterality);
            hitUsed[idx] = true;
         }
      // }    
      // std::cout << "BEST:\n";
      // for (auto& idx : cand.indices) {
      //    auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
      //    std::cout << hit->GetPlane() << '\t' << hit->GetLayer() << '\t' << hit->GetCell() << '\t' << hit->GetLaterality() << '\n';
      }
   }
}


void ConvDriftTubeRawData::FindLateralitySlope(const TClonesArray * hits, const std::vector<std::vector<int>>& clusters) {

   struct TrackResult {
      std::vector<int> latCombination {};
      double quality;
      std::vector<int> indices;
   };

   std::vector<TrackResult> candidates;

   const auto HCELL = DriftTubeDet->GetConfParF("DriftTube/cellHeight") / 10 + DriftTubeDet->GetConfParF("DriftTube/plateThickness") / 10; // Al plate thickness
   const auto WCELL = DriftTubeDet->GetConfParF("DriftTube/cellWidth") / 10; 
   
   int clusterID {-1};
   
   for (const auto& hitsIdx : clusters) {
      int nHits {hitsIdx.size()};
      if (nHits < 3) continue;
      
      TrackResult bestInCluster {{}, 1e9, hitsIdx};
      
      struct Point {
         double x;
         double y;
      };
      
      std::vector<Point> points;
      
      int nCombinations = 1 << 4;
      
      for (int i = 0; i != nCombinations; ++i) {
         std::vector<int> currentLats(4, 0); // if no hit in that layer, lat stays 0
         for (int idx : hitsIdx) {
            auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
            int layer {hit->GetLayer()};
            currentLats[layer] = ((i >> layer) & 1) ? 1 : -1;
            double y {(layer - 1.5) * HCELL};
            double x_wire {(hit->GetCell() + (layer % 2 == 1 ? 0.5 : 1.0)) * WCELL};
            double drift_space {(hit->GetTimestamp() - TPED) * VDRIFT}; // negative time hits??
            double x {x_wire + (currentLats[layer] * drift_space)};
            if (hit->GetTimestamp() < 0) {
               x = x_wire;
               currentLats[layer] = 0;
            } else if (hit->GetTimestamp() > WCELL / VDRIFT) {
               x = x_wire + currentLats[layer] * 0.5 * WCELL;
            }
            points.push_back({x, y});
         }
         
         std::vector<double> slopes {};
         double slopesMean {};
         double slopesStdev {};
         
         // std::cout << i << "\t size: " << points.size() << "\t hits: " << hitsIdx.size() << '\n';
         
         for (size_t j = 0; j != points.size(); ++j) {
            for (size_t k = j + 1; k < points.size(); ++k) {
               double slope {(points[j].y - points[k].y) / (points[j].x - points[k].x)};
               if (std::isnan(slope) || std::isinf(slope)) continue;
               // std::cout << i << "\t size: " << slopes.size() << "\t slope:" << slope << '\n';
               slopes.push_back(slope);
               slopesMean += slope;
            }
         }
         
         points.clear();
         
         slopesMean = slopesMean / slopes.size();
         
         // std::cout << i << "\t size: " << slopes.size() << "\t mean:" << slopesMean << '\n';
         
         for (const auto& slope : slopes) {
            slopesStdev += std::pow(slope - slopesMean, 2);
         }
         
         slopesStdev = slopesStdev / slopes.size();
         
         if (slopesStdev < bestInCluster.quality) {
            bestInCluster.latCombination = currentLats;
            bestInCluster.quality = slopesStdev;
            bestInCluster.indices = hitsIdx;
         }
         slopes.clear();
      }
      candidates.push_back(bestInCluster);
   }
   
   if (!candidates.empty()) {
      std::sort(candidates.begin(), candidates.end(), [](const TrackResult& a, const TrackResult& b) {return a.quality < b.quality;});
      
      std::vector<bool> hitUsed(static_cast<int>(hits->GetEntries()), false);
      
      for (auto& cand : candidates) {
         bool conflict {false};
         // if (cand.quality > 10) continue;
         ++clusterID;
         std::cout << "clusterID: " << clusterID << '\n';

         // for (int idx : cand.indices) {
         //    if (hitUsed[idx]) {
         //       break;
         //    } else {
         //       auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
         //       int layer {hit->GetLayer()};
   
         //       int lat {cand.latCombination[layer]};
   
         //       hit->setLaterality(lat);
         //       hit->setClusterID(clusterID);

         //       // std::cout << lat << "\t";
   
         //       hitUsed[idx] = true;
         //    }
         // }

         for (int idx : cand.indices) {
            if (hitUsed[idx]) {
               conflict = true;
               break;
            }
         }
   
         if (!conflict) {
            for (int idx : cand.indices) {
               auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
               int layer {hit->GetLayer()};
   
               int lat {cand.latCombination[layer]};
   
               hit->setLaterality(lat);
               hit->setClusterID(clusterID);

               std::cout << hit->GetClusterID() << '\t' << hit->GetPlane() << '\t' << hit->GetLayer() << '\t' << hit->GetCell() << '\t' << hit->GetLaterality() << '\n';
   
               hitUsed[idx] = true;
            }
         }

         // for (int idx : cand.indices) {
         //    auto hit = dynamic_cast<DriftTubeHit*>(hits->At(idx));
         //    std::cout << hit->GetClusterID() << '\t' << hit->GetPlane() << '\t' << hit->GetLayer() << '\t' << hit->GetCell() << '\t' << hit->GetLaterality() << '\n';
         // }
      }
   }

   
   candidates.clear();
}