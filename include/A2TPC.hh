/***** Time Projection Chamber
A Postuma 2021 *****/

#ifndef A2TPC_h
#define A2TPC_h 1

//necessary header files
#include "A2Target.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4Cache.hh"
#include "A2SD.hh"
#include "A2VisSD.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "A2ElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "CLHEP/Units/SystemOfUnits.h"


class A2TPC : public A2Target //define this class as part of class A2Target
{
public: //declare all public functions
	A2TPC(); //constructor
	~A2TPC(); //destructor

	virtual G4VPhysicalVolume* Construct(G4LogicalVolume *MotherLogical, G4double=0); //the volume it returns

	//functions called called by main constructor
	void MakeVessel();
	void MakeAnodeCathode();
	void MakeGrid();
	void MakeSensitiveDetector();
	void MakeField();
	void PlaceParts();
	void DefineMaterials();
	void ReadParameters(const char*);
	void SetIsOverlapVol(G4int isOv){ fIsOverlapVol = isOv; }

private: //private declarations

	//nist manager
	G4NistManager* fNistManager;
	
	//overlap checker
	G4int fIsOverlapVol;

	//sensitive detector stuff
	G4Region* fRegionAnode;
        A2SD* fAnodeSD;
        A2VisSD* fAnodeVisSD;

	//electric field stuff
	//G4ElectricField*        pEMfield;
        //G4FieldManager*         localFieldMgr;
        //G4EqMagElectricField*   pEquation;
        //G4ChordFinder*          pChordFinder;
	//G4UniformElectricField* fElectricField;
	//G4EqMagElectricField* fEquation;
	//G4VIntegrationDriver* fIntegrationDriver;
	//G4ChordFinder* fChordFinder;
	A2ElectricField* fElectricField;
	G4Region* fRegionActiveGas;

	//logical and physical volumes that are part of every detector class
	G4LogicalVolume* fMotherLogic; //logical volume of the mother
	G4LogicalVolume* fMyLogic;  //logical Volume for this detector
	G4VPhysicalVolume* fMyPhysi;   //physical volume for this detector
	
	//volumes specifically for this detector class
	//volumes in MakeVessel()
	G4LogicalVolume* fVesselLogic;
	G4LogicalVolume* fMainCellLogic;
	//G4LogicalVolume* fMainCellHeLogic;
	G4LogicalVolume* fMainCellEndLogic;
	G4LogicalVolume* fExtCellLogic;
	//G4LogicalVolume* fExtCellHeLogic;
	G4LogicalVolume* fConeLogic;
	G4LogicalVolume* fVesselHeLogic;
	//G4LogicalVolume* fConeHeLogic;
	G4LogicalVolume* fBeWindowLogic;
	//volumes in MakeAnodeCathode()
	G4LogicalVolume* fAnodeLogic;
	G4LogicalVolume* fAnodeCentreLogic;
        G4LogicalVolume* fAnodeRingLogic;
	G4Tubs* fAnodeSec[6]; //solid for each wire
        G4LogicalVolume* fAnodeSecLogic[6]; //logic of each wire

        //G4LogicalVolume* fAnodeSec1Logic;
        //G4LogicalVolume* fAnodeSec2Logic;
        //G4LogicalVolume* fAnodeSec3Logic;
        //G4LogicalVolume* fAnodeSec4Logic;
        G4LogicalVolume* fCathodeLogic;
	//volumes in MakeGrid()
	G4LogicalVolume* fGridLogic;
	//declare array longer than needed, as C++ can't process variable length arrays 
	G4Tubs* fWire[200]; //solid for each wire
        G4LogicalVolume* fWireLogic[200]; //logic of each wire

	//geometric parameters
	G4double fLength;
	G4double fRadius;
	G4double fThickness;
	G4double fConeLength;
	G4double fExtension;
	G4double fExtRadius;
	G4double fEndThickness;
	G4double fBeThickness;
	G4double fGThickness;
	G4double fAnodeRadius;
	G4double fAnodeDistance;
	G4int fRadialSecs;
	G4int fAngularSecs;
	G4double fAlThickness;
	G4double fCathodeRadius;
	G4double fCathodeDistance;
	G4double fWireThickness;
	G4double fWireSpacing;
	G4int fHeIsotope;
	G4double fHePressure;
	G4String fHeMaterial;

	//for placing anode
	//G4VPhysicalVolume** fAnodePhysi;
} ; 

#endif
