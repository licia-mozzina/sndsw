#include "DriftTubeHit.h"
#include "DriftTube.h"
#include "TROOT.h"
#include "FairRunSim.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include <TRandom.h>
#include <iomanip>

// -----   Default constructor   -------------------------------------------
DriftTubeHit::DriftTubeHit() 
   : TObject(),
    fDetectorID(-1)
{
}
// -----   Standard constructor   ------------------------------------------
DriftTubeHit::DriftTubeHit(Int_t detID) 
   : TObject(),
   fDetectorID(detID)
{
   flag = true;
}
// -----   Constructor from from MiniDT raw data  ------------------------------------------
// DriftTubeHit::DriftTubeHit(Int_t detID, const Double_t& m_timestamp, const Int_t &m_station, const Int_t &m_layer, const Int_t &m_wire) : fDetectorID(detID), timestamp(m_timestamp), station(m_station), layer(m_layer), wire(m_wire) 
// {
//    flag = true;
// }
DriftTubeHit::DriftTubeHit(Int_t detID, const Double_t& m_timestamp) : TObject(), fDetectorID(detID), timestamp(m_timestamp)
{
   flag = true;
}

// -----   constructor from point class  ------------------------------------------
DriftTubeHit::DriftTubeHit(int detID, std::vector<DriftTubePoint *> V, std::vector<Float_t> W)
{
   DriftTube *DriftTubeDet = dynamic_cast<DriftTube *>(gROOT->GetListOfGlobals()->FindObject("DriftTube"));
   // Float_t timeResol = DriftTubeDet->GetConfParF("DriftTube/timeResol"); // example

   // nSides = 1;
   // for (unsigned int j = 0; j < 16; ++j) {
   //    signals[j] = -1;
   //    times[j] = -1;
   // }

   for (auto p = std::begin(V); p != std::end(V); ++p) {

      Double_t signal = (*p)->GetEnergyLoss();
      // Find the distance from MCPoint to the center of cell (the anode)
      TVector3 vLeft, vRight;
      TVector3 impact((*p)->GetX(), (*p)->GetY(), (*p)->GetZ());
      DriftTubeDet->GetPosition(detID, vLeft, vRight);
      Double_t distance = (vLeft - impact).Perp(); // transverse component

      // for the timing - what to do?
      Double_t ptime = (*p)->GetTime();
   }
   // what needs to be set: distance? time?

   LOG(DEBUG) << "signal created";
}

// -----   Destructor   ----------------------------------------------------
DriftTubeHit::~DriftTubeHit() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void DriftTubeHit::Print()
{
   std::cout << "-I- DriftTubeHit: DriftTube hit " << " in station " << GetStation();
   if (isVertical()) {
      std::cout << " vertical plane ";
   } else {
      std::cout << " horizontal plane ";
   }
   std::cout << "layer nr " << GetLayer() << " cell nr " << GetCell() << std::endl;
}
// -------------------------------------------------------------------------

TVector3 DriftTubeHit::GetPosition() { // cambiare nome? fare altra funzione?
   TVector3 position {};
   std::string node = Form("/Detector_0/volDriftTubePlane_%d/volLayer_%d/volCell_%d/volAnode_2", GetPlane(), GetLayer(), GetDetectorID());
   auto navigator = gGeoManager->GetCurrentNavigator();

   if (gGeoManager->cd(node.c_str())) {
      TGeoShape *shape = gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
      // Check if the shape is a TGeoBBox (Box) and get dimensions
      if (shape->InheritsFrom("TGeoBBox")) {
         TGeoBBox *box = dynamic_cast<TGeoBBox *>(shape);
         const Double_t *origin = box->GetOrigin();  // gGeoManager->GetCurrentNode()->GetMatrix()->GetTranslation();

         navigator->cd(node.c_str());
         Double_t localPosition[3] {origin[0] + static_cast<double>((timestamp - TPED) * VDRIFT * laterality), origin[1], origin[2]};
         Double_t globalPosition[3] {};
         navigator->LocalToMaster(localPosition, globalPosition);

         position.SetXYZ(globalPosition[0], globalPosition[1], globalPosition[2]);

         // std::cout << node << '\n';
         // std::cout << "--- Geometry Debug for: " << node << " ---" << "\n";
         // std::cout << "local: " << localPosition[0] << "\t" << localPosition[1] << "\t" << localPosition[2] << "\n";
         // std::cout << "global: " << globalPosition[0] << "\t" << globalPosition[1] << "\t" << globalPosition[2] << "\n";
      }
   }

   return position;
}

ClassImp(DriftTubeHit)
