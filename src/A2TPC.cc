/***** Time Projection Chamber *****
 * Small TPC to be used in TPC+CB experiment.
 * Based on specifications by V Sokhoyan and E Maev.
 * Reads parameters from data/TPC.dat.
 * Requires class A2DriftModel for proper electron propagation.
 ***** AC Postuma 2021 *****/

#include "A2TPC.hh"

#include "G4SDManager.hh" //manage sensitive detector
#include "G4Tubs.hh" //cylinder
#include "G4Cons.hh" //cone
#include "G4NistManager.hh"  //element manager
#include "G4VisAttributes.hh" //visualization
#include "G4PVPlacement.hh" //placement
#include "G4UnionSolid.hh" //union of several solide
#include "G4RotationMatrix.hh" //rotate and place solids
#include "G4SystemOfUnits.hh" //units
#include "G4PhysicalConstants.hh" //constants
#include "A2SD.hh" //sensitive detector
#include "A2VisSD.hh" //visual sensitive detector
#include "A2ElectricField.hh" //my electric field class
#include "A2DriftModel.hh" //my electron propagation model
#include "G4FastSimulationManager.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

#include "G4FieldManager.hh" //manage the fields
#include "G4TransportationManager.hh" //transport through fields
#include "CLHEP/Units/SystemOfUnits.h" //units
#include "G4IntegrationDriver.hh"

using namespace CLHEP; //units

/***** Constuctor *****/
A2TPC::A2TPC(){
	//initiate pointers to NULL
	
	//first, volumes: the three that are in every target file
	fMotherLogic = NULL; //mother logical volume
	fMyLogic = NULL; //the constructed logical volume
	fMyPhysi = NULL; //the constructed physical volume
	
	//then the volumes specific to this target
	//volumes in MakeVessel()
	fVesselLogic=NULL; //contains all vessel parts
	fMainCellLogic = NULL; //main vessel shell
	fMainCellEndLogic = NULL; //flat end of main vessel
	fExtCellLogic = NULL; //extension cells
	fConeLogic = NULL; //conical end of main vessel
	fBeWindowLogic = NULL; //beryllium beam windows
	fVesselHeLogic = NULL; //helium inside vessel (union solid)
	//volumes in MakeAnodeCathode()
	fAnodeLogic = NULL; //contains all anode parts
        fAnodeCentreLogic= NULL; //centre of anode
        fAnodeRingLogic= NULL; //ring around anode centre
	fCathodeLogic= NULL; //cathode (one piece)
	for(G4int i=0;i<6;i++)fAnodeSecLogic[i]=NULL; //array of anode segments
	//volumes in MakeGrid()
	fGridLogic= NULL; //contains entire grid
	for(G4int j=0;j<200;j++)fWireLogic[j]=NULL; //array of grid wires
	
	//some sensitive detector stuff
        fAnodeSD = NULL; //sensitive detector
        fAnodeVisSD = NULL; //visualization
        fRegionAnode = new G4Region("Anode"); //set anode as a region

	//some E field stuff
	fElectricField = NULL; //electric field object
	fRegionActiveGas = new G4Region("ActiveGas"); //set helium as a region

	//initiate the nist manager
	fNistManager=G4NistManager::Instance();

	//initiate other relevant pointers here
	fIsOverlapVol=true; //check overlaps in geometry

	//read dimensions from a parameter file
	ReadParameters("data/TPC.dat");
}

/***** Destructor *****/
A2TPC::~A2TPC(){
	//remove things that requuire manual deletion
	if(fRegionAnode) delete fRegionAnode; //anode region
	if(fRegionActiveGas) delete fRegionActiveGas; //active gas region
	if(fAnodeSD) delete fAnodeSD; //sensitive detector
	if(fAnodeVisSD) delete fAnodeVisSD; //sensitive detector
	if(fElectricField) delete fElectricField; //electric field
}

/***** This function in called in DetectorSetup to build the TPC *****/
G4VPhysicalVolume* A2TPC::Construct(G4LogicalVolume* MotherLogical, G4double Z0){
	
	//get pointer to mother logical volume
	fMotherLogic=MotherLogical;
	
	//output a message
	G4cout<< "A2TPC::Construct() Building the TPC." <<G4endl;
	
	//call functions
	DefineMaterials(); //define relevant materials
	MakeVessel(); //build fVesselLogic
	MakeAnodeCathode(); //build fAnodeLogic and fCathodeLogic
	MakeGrid(); //build fGridLogic
	MakeSensitiveDetector(); //assign anode to sensitive detector
	MakeField(); //create electric field inside the active volume
	PlaceParts(); //place the seperate detector parts into fMyLogic
	
	//place fMyLogic into fMotherLogic
	fMyPhysi = new G4PVPlacement(0, G4ThreeVector(0,0,Z0) ,fMyLogic,"TPC",fMotherLogic,false,1,fIsOverlapVol);

	//and return this physical volume describing the target!
	return fMyPhysi;
}

/***** This function builds the target vessel and fills the target with helium *****/
void A2TPC::MakeVessel(){
	/***** shapes ******/
	//overall vessel shape
	G4Tubs *fVessel = new G4Tubs("Vessel", //name
				0, //inner radius
				fRadius+fThickness, //outer radius
				fLength/2+fExtension+fConeLength, //half length
				0.*deg, //start angle
				360.*deg); //spanning angle
	//main cell steel shell
	G4Tubs *fMainCell = new G4Tubs("MainCell",
				fRadius, 
				fRadius+fThickness,
				fLength/2,
                                0.*deg,
                                360.*deg);
	//circular end of main cell, with hole for extension
	G4Tubs *fMainCellEnd = new G4Tubs("MainCellEnd",
				fExtRadius+fThickness,
				fRadius+fThickness,
				fEndThickness/2,
				0.*deg,
				360.*deg);
	//conical main cell end on the other side
        G4Cons *fMainCellCone = new G4Cons("MainCellCone",
                                fRadius, //inner radius at one end
                                fRadius+fThickness, //outer
                                fExtRadius, //inner radius at other end
                                fExtRadius+fThickness, //outer
                                fConeLength/2, //half length
                                0.*deg, //start angle
                                360.*deg); //spanning angle
	//helium inside main cell
	G4Tubs *fMainCellHe = new G4Tubs("MainCellHelium",
				0,
				fRadius,
				fLength/2,
				0.*deg,
				360.*deg);
        //helium inside conical part of cell
	G4Cons *fConeHe = new G4Cons("ConeHe",
                                0,
                                fRadius,
                                0,
                                fExtRadius,
                                fConeLength/2,
                                0.*deg,
                                360.*deg);
	//extension tube
	//made once and placed twice
	G4Tubs *fExtCell = new G4Tubs("ExtCell",
				fExtRadius,
				fExtRadius+fThickness,
				fExtension/2,
				0.*deg,
				360.*deg);
	//helium inside main extension tubes
	G4Tubs *fExtCellHe1 = new G4Tubs("ExtCellHe",
				0,
				fExtRadius,
				fExtension/2,
				0.*deg,
				360.*deg);
	G4Tubs *fExtCellHe2 = new G4Tubs("ExtCellHe",
				0,
				fExtRadius,
				fExtension/2+fEndThickness/2,
				0.*deg,
				360.*deg);
	//beryllium window at the end of each extension tube
        G4Tubs *fBeWindow = new G4Tubs("BeWindow",
                                0,
                                41.*mm,
                                fBeThickness/2,
                                0.*deg,
                                360.*deg);
	
	/***** create union solid of all helium *****/	
	G4UnionSolid *fHeUnion1 = new G4UnionSolid("HeUnion1",
			fMainCellHe,
			fConeHe,
			0,
			G4ThreeVector(0,0,(fLength+fConeLength)/2));
	G4UnionSolid *fHeUnion2 = new G4UnionSolid("HeUnion2",
			fHeUnion1,
			fExtCellHe1,
			0,
			G4ThreeVector(0,0,(fLength+fExtension)/2+fConeLength));
	G4UnionSolid *fVesselHe = new G4UnionSolid("VesselHe",
			fHeUnion2,
			fExtCellHe2,
			0,
			G4ThreeVector(0,0,-(fLength+fExtension+fEndThickness)/2));
	
	/**** logical volumes ******/
	//vessel
	fVesselLogic = new G4LogicalVolume
		(fVessel,
		 fNistManager->FindOrBuildMaterial("G4_AIR"),
		 "VesselLogic");
	//main cell
	fMainCellLogic = new G4LogicalVolume
            	(fMainCell,   //solid
            	 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "MainCellLogic"); //name for logical volume
	//main cell endpiece
	fMainCellEndLogic = new G4LogicalVolume
		(fMainCellEnd,
		 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
		 "MainCellEndLogic");
	//conical end to main cell
        fConeLogic = new G4LogicalVolume
                (fMainCellCone,
                 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
                 "ConeLogic");
	//helium
	fVesselHeLogic = new G4LogicalVolume
		(fVesselHe,
		 fNistManager->FindOrBuildMaterial(fHeMaterial),
		 "VesselHeLogic");
	//extension cells
	fExtCellLogic = new G4LogicalVolume
		(fExtCell,
		 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
		 "ExtCellLogic");
	//beryllium windows
	fBeWindowLogic = new G4LogicalVolume
                (fBeWindow,
                 fNistManager->FindOrBuildMaterial("ATBerylliumW"),
                 "BeWindowLogic");
	

	/**** set visualization attributes *****/
	G4VisAttributes* lblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
    	G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
    	G4VisAttributes* cyan   = new G4VisAttributes( G4Colour(0.0,1.0,1.0,0.5)  );
	fVesselLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
	fMainCellLogic->SetVisAttributes(grey);
	//fMainCellLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fVesselHeLogic->SetVisAttributes(cyan);
	//fVesselHeLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fMainCellEndLogic->SetVisAttributes(grey);
    	fConeLogic->SetVisAttributes(grey);
	//fMainCellEndLogic->SetVisAttributes(G4VisAttributes::Invisible);
    	//fConeLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fExtCellLogic->SetVisAttributes(grey);
	//fExtCellLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fBeWindowLogic->SetVisAttributes(lblue);
	//fBeWindowLogic->SetVisAttributes(G4VisAttributes::Invisible);

	/**** place things inside of VesselLogic ****/
	//main cell
    	new G4PVPlacement (0, //rotation
                     G4ThreeVector(0,0,0), //position
                     fMainCellLogic, //logical volume
                     "MainCellPlacement", //name
                     fVesselLogic, //mother logic
                     false, //pMany = false always
                     2,fIsOverlapVol); //unique copy number
	//main cell endpiece
    	new G4PVPlacement(0,
		    G4ThreeVector(0,0,-(fLength+fEndThickness)/2),
		    fMainCellEndLogic,
		    "MainCellEndPlacement",
		    fVesselLogic,
		    3,fIsOverlapVol);
    	//main cell cone
	new G4PVPlacement(0,
		    G4ThreeVector(0,0,(fLength+fConeLength)/2),
		    fConeLogic,
		    "ConePlacement",
		    fVesselLogic,
		    4,fIsOverlapVol);
	//extension cells
	new G4PVPlacement(0,
		   G4ThreeVector(0,0,(fLength+fExtension)/2+fConeLength),
		   fExtCellLogic,
	  	   "ExtCellPlacement1",
		   fVesselLogic,
		   6,fIsOverlapVol);
	new G4PVPlacement(0,
		   G4ThreeVector(0,0,-(fLength+fExtension)/2-fEndThickness),
	 	   fExtCellLogic,
		   "ExtCellPlacement2",
		   fVesselLogic,
		   7,fIsOverlapVol);
	//beryllium windows - for now overlapping with helium
        new G4PVPlacement(0,
                   G4ThreeVector(0,0,fLength/2 + fExtension -fBeThickness/2+fConeLength),
                   fBeWindowLogic,
                   "BeWindowPlacement1",
                   fVesselLogic,
                   10,fIsOverlapVol);
        new G4PVPlacement(0,
                   G4ThreeVector(0,0,-(fLength/2+fExtension-fBeThickness/2+fEndThickness)),
                   fBeWindowLogic,
                   "BeWindowPlacement2",
                   fVesselLogic,
                   11,fIsOverlapVol);
}

/***** This function creates a segmented anode and a cathode *****/
void A2TPC::MakeAnodeCathode(){
        /***** solids for anode geometry *****/
        //main volume to hold sub-pieces
	G4Tubs *fAnode = new G4Tubs("Anode", //name
                                        0, //inner radius
                                        fRadius, //outer radius
                                        (fGThickness)/2, //half length
                                        0.*deg, //start angle
                                        360.*deg); //spanning angle
	
	//circular central piece (G-10)
        G4Tubs *fAnodeCentre = new G4Tubs("AnodeCentre",
                                        0,
                                        fRadPad*mm,
                                        fGThickness/2,
                                        0.*deg,
                                        360.*deg);
	//ring around central piece (G-10)
        G4Tubs *fAnodeRing = new G4Tubs("AnodeRing", 
					fRadPad*mm,
                                        fRadRing*mm,
					fGThickness/2,
                                        0.*deg,
                                        360.*deg);
	
	/***** logical volumes *****/	
	//volume holding entire anode (helium)
	fAnodeLogic = new G4LogicalVolume
                        (fAnode,
                         fNistManager->FindOrBuildMaterial(fHeMaterial),
                         "AnodeLogic");
	
	//circular central piece (G-10)
        fAnodeCentreLogic = new G4LogicalVolume
                        (fAnodeCentre,
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeCentreLogic");
	//ring around central piece (G-10)
        fAnodeRingLogic = new G4LogicalVolume
                        (fAnodeRing,
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeRingLogic");
	
	/***** visualization attributes *****/
        G4VisAttributes* lblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
        G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
        fAnodeLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
        fAnodeCentreLogic->SetVisAttributes(lblue);
        fAnodeRingLogic->SetVisAttributes(lblue);

        /***** place anode parts to fAnodeLogic *****/
        //circular central piece (G-10)
	new G4PVPlacement(0, //rotation
                        G4ThreeVector(0,0,0), //placement
                        fAnodeCentreLogic, //logical volume
                        "AnodeCentrePlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        //1);
			2+4*(fAngularSecs),fIsOverlapVol); //unique copy number
        //ring around central piece (G-10)
	new G4PVPlacement(0, //rotation
                        G4ThreeVector(0,0,0), //placement
                        fAnodeRingLogic, //logical volume
                        "AnodeRingPlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        //2);
			1+4*(fAngularSecs),fIsOverlapVol); //unique copy number
	
	/*** segments of anode: reads specifications from parameter file ****/
	/***** define solid and logical volume for each radial section *****/
	//new design: add some 
	G4double radii[5]={fRadRing,fRad1,fRad2,fRad3,fRad4};
	for(G4int k=0; k<4; k++){
		//section of each row (G-10)
		fAnodeSec[k] = new G4Tubs("AnodeSec",
                                        radii[k],
					radii[k+1],
					fGThickness/2,
                                        0.*deg,
					360.*deg/fAngularSecs);
        	fAnodeSecLogic[k] = new G4LogicalVolume
                        (fAnodeSec[k],
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeSecLogic");
		fAnodeSecLogic[k]->SetVisAttributes(lblue);

		/***** place ring segments to anode *****/
        	for(G4int h=0; h<fAngularSecs;h++){ //for however many angular sections were defined in paramter file
          		//define the angle
                	G4double theta = h*360.*deg/fAngularSecs;
                	//generate a rotation matrix
                	G4RotationMatrix* rot = new G4RotationMatrix(theta,0,0);
                	//place one section of each ring at this angle
			//section of first row (G-10)
                	new G4PVPlacement(rot, //rotation
                        G4ThreeVector(0,0,0),
                        fAnodeSecLogic[k], //logical volume
                        "AnodeSecPlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        1+h*4+k,fIsOverlapVol);
			//1+h+k*fAngularSecs,fIsOverlapVol); //unique copy number
		}
	}
        
	/***** build solids for cathode geometry *****/
        //aluminum cathode
	G4Tubs *fCathode = new G4Tubs("Cathode",
                                        0,
                                        fCathodeRadius,
                                        fAlThickness/2,
                                        0.*deg,
                                        360.*deg);
        
	/***** cathode logical volumes *****/
        //aluminum cathode
	fCathodeLogic = new G4LogicalVolume(fCathode,
                        fNistManager->FindOrBuildMaterial("Aluminum"),
                        "CathodeLogic");
        
	/***** visulization attributes *****/
	fCathodeLogic->SetVisAttributes(lblue);
}

/**** This function builds the wire grid *****/
void A2TPC::MakeGrid(){
	/***** make main volume containing wires *****/
	//solid: cylinder
	G4Tubs* fGrid = new G4Tubs("Grid", //name
				0, //inner radius
				fRadius, //outer radius
				2.*mm, //half length in Z
				0, //start angle
				360.*deg); //spanning angle
	//logical volume
	fGridLogic = new G4LogicalVolume(
			fGrid, //solid
			fNistManager->FindOrBuildMaterial(fHeMaterial), //material: helium, same as fVesselHeLogic
			"GridLogic"); //name

	/***** define parameters to make wires ******/
	G4int nWires = fRadius/fWireSpacing-1; //how many different length wires to build
	G4double fHWL[200]; //half length of each wire
	G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
        fGridLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
	
	/***** make and place each individual wires ******/
	for(G4int l=0; l<nWires; l++){
		//find length needed to cross cylinder at a certain dinstace from the origin
		fHWL[l] = sqrt(fRadius*fRadius - l*l*fWireSpacing*fWireSpacing);
		//define wire with correct length
		fWire[l] = new G4Tubs("Wire",
			0, //inner radius
			fWireThickness, //outer radius
			fHWL[l], //half length
			0.*deg, //start angle
			360.*deg); //spanning angle
		//make the logical volume
		fWireLogic[l] = new G4LogicalVolume(fWire[l],
			fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
			"WireLogic");
        	fWireLogic[l]->SetVisAttributes(grey);
		//determine distance from origin
		G4double XYPos = l*fWireSpacing;
		//create rotation matrices to align wire with X or Y axis
		G4RotationMatrix *xrot = new G4RotationMatrix(0,90.*deg,0);
		G4RotationMatrix *yrot = new G4RotationMatrix(90.*deg,90.*deg,0);
		//place 4 copies of each wire:
		//pointing along Y, displaced from origin in +X
		new G4PVPlacement(
				xrot, //rotation
				G4ThreeVector(XYPos,0,fWireThickness), //spatial coordinates 
				fWireLogic[l], //logical volume
				"WirePlacementX", //name
				fGridLogic, //mother logical
				false, //pMany: always false
				143, //unique copy no
			       	false); //check overlaps: not here
		//pointing along Y, displaced from origin in +X
		new G4PVPlacement(
				xrot,
				G4ThreeVector(-XYPos,0,fWireThickness),
				fWireLogic[l],
				"WirePlacementX2",
				fGridLogic,
				false,
				143, false);
		//pointing along X, displaced from origin in +Y
		new G4PVPlacement(
				yrot,
				G4ThreeVector(0,XYPos,-fWireThickness),
				fWireLogic[l],
				"WirePlacementY",
				fGridLogic,
				false,
				144,false);
		//pointing along X, displaced from origin in -Y
		new G4PVPlacement(
				yrot,
				G4ThreeVector(0,-XYPos,-fWireThickness),
				fWireLogic[l],
				"WirePlacementY2",
				fGridLogic,
				false,
				144,false);
	}

}

/***** this function makes the anode into a sensitive detector *****/
void A2TPC::MakeSensitiveDetector(){
	/***** define and register sensitive detector *****/
        if(!fAnodeSD){ //if SD not already defined
        G4SDManager* SDman = G4SDManager::GetSDMpointer(); //get pointer to SD manager
        fAnodeSD = new A2SD("AnodeSD",2+4*fAngularSecs); //create a new SD
        SDman->AddNewDetector(fAnodeSD); //add this detector to the SD manager

        /***** set each piece of anode as part of the sensitive detector *****/
	
	for(G4int n=0; n<4; n++){ //for each anode section
		fAnodeSecLogic[n]->SetSensitiveDetector(fAnodeSD); //make sensitive
		fRegionAnode->AddRootLogicalVolume(fAnodeSecLogic[n]); //add to region
	}
	fAnodeCentreLogic->SetSensitiveDetector(fAnodeSD); //make sensitive
        fAnodeRingLogic->SetSensitiveDetector(fAnodeSD); //make sensitive
        fRegionAnode->AddRootLogicalVolume(fAnodeCentreLogic); //add to region
	fRegionAnode->AddRootLogicalVolume(fAnodeRingLogic); //add to region
    }
}

/***** this function creates electric field and field processes for the TPC *****/
void A2TPC::MakeField(){
	   /***** create uniform electric field in helium ****/
	   // fElectricField = new A2ElectricField(); //create (zero) electric field
	   // f qqfElectricField->Construct(2); //set the field strength in kV/cm

	/***** attach the active gas region without Geant4 field transport ****/
	// The drift model below parametrizes electron transport. Installing a
	// Geant4 electric field here also steers all charged particles in the TPC
	// gas and can make event tracking effectively stall for Compton input.

	fRegionActiveGas->AddRootLogicalVolume(fVesselHeLogic); //create active gas region

	/***** set specific production cuts for active gas region *****/
	G4ProductionCuts* TPCcuts = new G4ProductionCuts(); //create custom cut
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(7*eV,1*GeV); //optimized at 7 eV lower limit
	//this gives the most linear charge-energy relationship
	fRegionActiveGas->SetProductionCuts(TPCcuts); //assign this cut to this region
	
	/***** attach model of electron drift to active gas region *****/
	A2DriftModel *driftPhys = new A2DriftModel("Electron Drift Model",fRegionActiveGas,this,fAnodeSD); //create model
	G4FastSimulationManager* driftMan = new G4FastSimulationManager(fRegionActiveGas); //call fast simulation manager
	driftMan->AddFastSimulationModel(driftPhys); //add model as a fast simulation
	driftMan->ActivateFastSimulationModel("driftPhys"); //activate model
}


/***** this function places all parts of the target to fMyLogic *****/
void A2TPC::PlaceParts(){
	/***** define fMyLogic *****/
	//create a simple cylinder with adequate room for main and extension cells
	G4Tubs *fMyLogicCyl = new G4Tubs
            ("MyLogicCyl",
             0,             //inner radius
             fRadius+fThickness,   //outer radius
             fLength/2+fExtension+25.*mm,  //length/2
             0.*deg,        //start angle
             360.*deg);     //final angle
	//turn it into a logical volume
	fMyLogic = new G4LogicalVolume
            (fMyLogicCyl,  //solid
             fNistManager->FindOrBuildMaterial("G4_AIR"), //material
             "MyLogic");    //name
	//and set it invisible
    	fMyLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
    
	/**** place anode and cathode to fVesselHeLogic *****/ 
        //first anode, close to beam entrance
	new G4PVPlacement(0,
                        G4ThreeVector(0,0,-fLength/2+fAnodeDistance),
                        fAnodeLogic,
                        "ANODE",
                        fVesselHeLogic,
                        30,fIsOverlapVol);
	//then cathode, close to beam exit
        new G4PVPlacement(0,
                        G4ThreeVector(0,0,fLength/2-fCathodeDistance),
                        fCathodeLogic,
                        "CATHODE",
                        fVesselHeLogic,
                        40,fIsOverlapVol); 
	/**** place grid to fVesselHeLogic ****/ 
	new G4PVPlacement(0,
			G4ThreeVector(0,0,-fLength/2+fAnodeDistance+3.*mm),
			fGridLogic,
			"GRID",
			fVesselHeLogic,
			50,fIsOverlapVol);
	/**** place helium containing anode, cathode, and grid to fVesselLogic ****/
	new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			fVesselHeLogic,
			"HELIUM",
			fVesselLogic,
			9,fIsOverlapVol);
	/**** place vessel inside fMyLogic *****/
	new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			fVesselLogic,
			"VESSEL",
			fMyLogic,
			10,fIsOverlapVol);

}

/***** this function reads target dimensions from a data file *****/
void A2TPC::ReadParameters(const char* file){
	//define keys contained in file
        char* keylist[] = { (char*) "TPC-Dim:", (char*) "Anode-Dim:", (char*) "Cathode-Dim:", (char*) "Grid-Dim:", (char*) "Helium:", (char*) "Run-Mode:", NULL};
        enum { ETPC_dim, EAnode_dim, ECathode_dim, EGrid_dim, EHelium, ERun_Mode, ENULL };
        //define variables needed for reading the file
        G4int ikey, iread;
        G4int ierr = 0;
        char line[256];
        char delim[64];
        //open the file
        FILE* pdata;
        if( (pdata = fopen(file,"r")) == NULL ){
                printf("Error opening detector parameter file: %s\n",file);
                return;
        }
        //start reading the file
        while( fgets(line,256,pdata) ){
        if( line[0] == '#' ) continue; //skip commented lines
        printf("%s\n",line);
        sscanf(line,"%s",delim);
        for(ikey=0; ikey<ENULL; ikey++) //look for the defined keys
            if(!strcmp(keylist[ikey],delim)) break;
        switch( ikey ){
        default:
	    //stop running if the file contains something strange
            printf("Unrecognised delimiter: %s\n",delim);
	    ierr++;
            break;
        case ETPC_dim: //dimensions of main target vessel
            iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf", //format of line
                           &fLength,&fRadius,&fThickness,&fConeLength, //pointers to associated variables
                           &fExtension,&fExtRadius,&fEndThickness,&fBeThickness);
            if( iread != 8) ierr++;
            break;
        case EAnode_dim: //dimensions of anode
            iread = sscanf(line,"%*s%lf%lf%i%lf%lf%lf%lf%lf%lf",
                            &fGThickness,&fAnodeDistance,&fAngularSecs,&fRadPad,&fRadRing,&fRad1,&fRad2,&fRad3,&fRad4);
            if (iread !=9) ierr++;
            break;
        case EGrid_dim: //dimensions of anode
            iread = sscanf(line,"%*s%lf%lf",
                            &fWireThickness,&fWireSpacing);
            if (iread !=2) ierr++;
            break;
        case ECathode_dim: //dimensions of cathode
            iread = sscanf(line,"%*s%lf%lf%lf",
                            &fAlThickness,&fCathodeRadius,&fCathodeDistance);
            if (iread !=3) ierr++;
            break;
	case EHelium:
	    iread = sscanf(line,"%*s%i%lf",
			    &fHeIsotope,&fHePressure);
	    if (iread !=2) ierr++;
	    break;
        case ERun_Mode: //run mode
            iread = sscanf(line,"%*s%d",
                           &fIsOverlapVol);
            if( iread != 1 ) ierr++;
            break;
	    }
        if( ierr ){ //if something strange is in the file, error and exit
    		printf("Fatal Error: invalid read of parameter line %s\n %s\n",
                   keylist[ikey],line);
            exit(-1);
                }
        }
}


/***** this function defines the materials used to build the target ******/
void A2TPC::DefineMaterials()
{
	G4double density, fractionmass;
	G4int ncomponents;
	//G4double pressure, temperature, a, z;

	/***** beryllium used in target windows *****/
	G4Material* BerylliumW = new G4Material("ATBerylliumW", density = 1.8480*g/cm3, ncomponents = 7);
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(14), 0.06*perCent);      //Si
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(4), 98.73*perCent);      //Be
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(6), 0.15*perCent);       //C
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(8), 0.75*perCent);       //O
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(12), 0.08*perCent);      //Mg
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(13), 0.1*perCent);       //Al
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(26), 0.13*perCent);      //Fe

	/***** anode and cathode materials ****/
        //copied from online source
        //Anode G-10
        G4Material* G10 = new G4Material("G-10", density= 1.700*g/cm3, ncomponents=4);
        G10->AddElement(fNistManager->FindOrBuildElement(14), 1); //silicon
        G10->AddElement(fNistManager->FindOrBuildElement(8), 2); //oxygen
        G10->AddElement(fNistManager->FindOrBuildElement(6), 3); //carbon
        G10->AddElement(fNistManager->FindOrBuildElement(1), 3); //hydrogen
        //Anode copper
        //G4Material* Cu = new G4Material("Copper",density= 8.960*g/cm3,ncomponents=1);
        //Cu->AddElement(fNistManager->FindOrBuildElement(29),1);
        //Cathode aluminum
        G4Material* Al = new G4Material("Aluminum", density= 2.700*g/cm3,ncomponents=1);
        Al->AddElement(fNistManager->FindOrBuildElement(13),1);
    
	
	/***** Active gas from A2ActiveHe3.cc *****/
	//Gas mixture and He3 management----------------------------------------------------

	//defining He3 according to
	// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/hadronprocess/731/1.html
	G4Element* ATHe3 = new G4Element("ATHe3", "ATHe3", ncomponents=1);
	ATHe3->AddIsotope((G4Isotope*)fNistManager->FindOrBuildElement(2)->GetIsotope(0), //isot. 0 = 3he
                      100.*perCent);
	G4Element* ATH = new G4Element("ATH","ATH",ncomponents=1);
	ATH->AddIsotope((G4Isotope*)fNistManager->FindOrBuildElement(1)->GetIsotope(0),100.*perCent);

	//See HeGasDensity.nb. Might be an idea to include N2 effect to density
	//G4double he3density = 0.00247621*g/cm3; //20bar, calculated from ideal gas law
	//G4double he3density = 0.0033*g/cm3; //20bar, calculated from ideal gas law
	//G4double he3density = 0.004125*g/cm3; //25bar, calculated from ideal gas law
	//G4double he3density = 0.00495*g/cm3; //30bar, calculated from ideal gas law
	G4double he3density = (fHePressure/20)*0.0033*g/cm3; //scale Phil's IG calculation according to pressure from parameter file 

	G4Material* GasMix = new G4Material("ATGasMix", he3density, ncomponents = 2,kStateGas,CLHEP::STP_Pressure,fHePressure*bar);
	GasMix->AddElement(ATHe3, 99.95*perCent);                                       //He3
	// GasMix->AddElement(fNistManager->FindOrBuildElement(2), 99.95*perCent);         //He4
	GasMix->AddElement(fNistManager->FindOrBuildElement(7), 0.05*perCent);           //N
	//Gas mixture and He3 end-----------------------------------------------------
	//----! IMPORTANT! Epoxy CURRENTLY TAKEN FROM A2 SIMULATION,------------------
	//NO IDEA WHETHER IT IS CORRECT OR NOT

	G4Material* He3GasPure = new G4Material("He3GasPure", he3density, ncomponents = 1,kStateGas,CLHEP::STP_Temperature,fHePressure*bar);
	He3GasPure->AddElement(ATHe3, 100.*perCent);

	//Active gas - has 10% hydrogen
	//density decreases accordingly: 0.9*1 + 0.1 * (2/3) = 0.96 of original density
	G4Material* He3ActiveGas = new G4Material("He3ActiveGas",he3density*0.96, ncomponents =2, kStateGas, CLHEP::STP_Temperature,fHePressure*bar);
	He3ActiveGas->AddElement(ATHe3,90.*perCent);
	He3ActiveGas->AddElement(ATH,10.*perCent);

	//define He4 in a similar manner
	G4Element* ATHe4 = new G4Element("ATHe4","ATHe4",ncomponents=1);
	ATHe4->AddIsotope((G4Isotope*)fNistManager->FindOrBuildElement(2)->GetIsotope(1),100.*perCent);
	G4double he4density = (fHePressure/20)*0.0033*g/cm3*(4/3); //scale Phil's IG calculation according to pressure from parameter file
	//4He edit: assume same number density as helium-3, but more nucleons means more mass

	G4Material* He4GasPure = new G4Material("He4GasPure",he4density, ncomponents=1,kStateGas,CLHEP::STP_Temperature,fHePressure*bar);
	He4GasPure->AddElement(ATHe4, 100.*perCent);


	G4Material* He4ActiveGas = new G4Material("He4ActiveGas",he4density*0.96, ncomponents =2, kStateGas, CLHEP::STP_Temperature,fHePressure*bar);
	He4ActiveGas->AddElement(ATHe4,90.*perCent);
	He4ActiveGas->AddElement(ATH,10.*perCent);
	//decide which version of helium you are using
	//G4String fHeMaterial;
	if(fHeIsotope==3){fHeMaterial="He3ActiveGas";
	} else if(fHeIsotope==4){fHeMaterial="He4ActiveGas";
	} else {fHeMaterial="ATGasMix";} //make the mix for any non-3,4 argument
}
