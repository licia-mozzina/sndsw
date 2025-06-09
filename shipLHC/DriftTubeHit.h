#ifndef SHIPLHC_DRIFTTUBEHIT_H_
#define SHIPLHC_DRIFTTUBEHIT_H_

#include "SndlhcHit.h"
#include "DriftTubePoint.h"
#include "TObject.h"
#include "TVector3.h"

class DriftTubeHit : public SndlhcHit {
public:
   /** Default constructor **/
   DriftTubeHit();
   /** Copy constructor **/
   DriftTubeHit(const DriftTubeHit &hit) = default;
   DriftTubeHit &operator=(const DriftTubeHit &hit) = default;
   explicit DriftTubeHit(Int_t detID);
   //  Constructor from MiniDT raw data
   explicit DriftTubeHit(Int_t detID, const Double_t& m_timestamp, const Int_t &m_station, const Int_t &m_layer, const Int_t &m_wire);
   //  Constructor from DriftTubePoint
   DriftTubeHit(int detID, std::vector<DriftTubePoint *>, std::vector<Float_t>);

   /** Destructor **/
   virtual ~DriftTubeHit();

   /** Output to screen **/
   void Print();

   // Overwrite functions from the parent SndlhcHit class
   // that don't make sense for the drift tubes
   Float_t GetSignal(Int_t nChannel = 0) { return -1; }
   Int_t GetnSiPMs() const { return -1; }

   void setInvalid() { flag = false; }
   bool isValid() const { return flag; }
   Double_t GetTimestamp() {return timestamp; }
   Int_t GetSystem() { return floor(fDetectorID / 10000); }
   Int_t GetStation() { return 1; } // do we need such a method? FIXME
   // Int_t GetStation() { return m_station; }
   // Int_t GetPlane() { return int(fDetectorID / 1000) % 10; }
   Int_t GetPlane() { return station; }
   // Int_t GetLayer() { return int(fDetectorID % 1000) / 100; }
   Int_t GetLayer() { return layer; }
   // Int_t GetCell() { return int(fDetectorID % 100); }
   Int_t GetCell() { return wire; }
   bool isVertical() { return GetPlane() == 1; } 
   // Double_t GetDistance(); // FIXME format of output

private:
   Float_t flag; ///< flag

   Double_t timestamp; // questa si droppa dopo matching? dopo ci sarà pedestal (variabile) con SND, tanto vale tenerla
   Int_t station, layer, wire;

   ClassDef(DriftTubeHit, 1);
};

#endif // SHIPLHC_DRIFTTUBEHIT_H_
