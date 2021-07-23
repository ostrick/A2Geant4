/***** Electron drift model for TPC *****
 * Implementation of TPC physics simulated in Heed/Garfield/Degrad.
 * Based off method described in https://arxiv.org/pdf/1806.05880.pdf and
 * https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE.
 * Uses G4VFastSimulation to implement manual definition of electron drift.
 ***** AC Postuma 2021 *****/

#include "A2HeedModel.hh"

#include "G4VPhysicalVolume.hh" //volumes
#include "G4Electron.hh" //particle the model applies to
#include "G4Gamma.hh" //another relevant (?) particle
#include "G4SystemOfUnits.hh" //units
#include "A2DetectorConstruction.hh" //detector construction
#include "G4RunManager.hh" //run
//#include <stdio.h>
#include "A2Target.hh" //targets
#include "A2SD.hh" //sensitive detector
//#include "DriftLineTrajectory.hh"
#include "G4TrackingManager.hh" //tracking of particles
#include "G4EventManager.hh" //events
#include "G4TransportationManager.hh" //particle transport
#include "G4VVisManager.hh" //track visualization

/**** Constructor *****/
A2HeedModel::A2HeedModel(G4String modelName, G4Region* actVol, A2Target* target, A2SD* anode)
: G4VFastSimulationModel(modelName, actVol), fA2Target(target), fA2SD(anode) { //fast simulation implements user defined physics response
	//initiate pointers
	fFakeStep = new G4Step(); //step used to call hit in SD
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint(); //step point
  	fFakePostStepPoint = fFakeStep->GetPostStepPoint(); //step point
  	fTouchableHandle   = new G4TouchableHistory(); //touchable for step
  	fpNavigator        = new G4Navigator(); //navigator to find SD
	fNaviSetup = false; //if setup has already been done
}

/***** Destructor *****/
A2HeedModel::~A2HeedModel(){
	//remove objects requiring manual deletion
	delete fFakeStep;
  	delete fpNavigator;
}

/***** Called in SteppingAction: checks particle type if model is applicable ******/
G4bool A2HeedModel::IsApplicable(const G4ParticleDefinition& particleType){
	G4String particleName = particleType.GetParticleName();
	if(particleName=="e-")return true; //only applicable to electrons
	return false;
}

/***** Called in SteppingAction: conditions in which to trigger model *****/
G4bool A2HeedModel::ModelTrigger(const G4FastTrack& fastTrack){
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	if (ekin <=1)return true; //trigger for kinetic energy below 1 keV
	return false;
}

/***** This function contains the main operation of the model *****/
void A2HeedModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep){
	//get all relevant step data from the track
	G4ThreeVector direction = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
	G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition()/cm;
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
	G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	//pass step data to transportation 
	Transport(fastStep, fastTrack, particleName, ekin, time, worldPosition.x(), worldPosition.y(), worldPosition.z(), direction.x(), direction.y(), direction.z());
}

//take just the delta electron run instructions
void A2HeedModel::Transport(G4FastStep& fastStep,const G4FastTrack& fastTrack, G4String particleName, double ekin_keV, double t, double x_cm, double y_cm, double z_cm, double dx, double dy, double dz){
	//G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	
	//calculate track parameters
	//everything here is incorrect - I'm just trying to have the simulation do something, then I'll get it to do the right thing
	//idea: first compute drift velocity. then drift the particle until a) it hits anode, or b) it hits the walls
	//find time taken to drift
	//x,y positions will be initial velocity times travel time, although some energy will be lost to the gas
	//z position will be drift velocity times travel time
	//might need to do a loop and keep processing it?
	G4double drift_vel = 200000; //cm/s... this is a lot... citing this paper here: https://journals.aps.org/pr/pdf/10.1103/PhysRev.117.470
	//probably need to implement a for loop here, or some alternate kind of stepping action
	//consider sideways motion to be negligible as a first approximation and just drift electrons to the anode	
	G4double x_pos = x_cm*mm;
	G4double y_pos = y_cm*mm;

	G4double pathLength = -115.5 - z_cm; //propogate to anode Z position
	//G4double radius = sqrt(x_cm*x_cm + y_cm+y_cm); //calculate radius - will it hit the anode?
	G4double time = pathLength/drift_vel; //time taken to drift to stopping place


	
	G4double z_pos = -115.5*mm;
	//G4double z_pos = -113*mm;
	G4ThreeVector position = G4ThreeVector(x_pos,y_pos,z_pos);

	//create any necessary secondary particles
	//fastStep.CreateSecondaryTrack();

	//set final track parameters
	//fastStep.KillPrimaryTrack(); //kill the step
	fastStep.SetPrimaryTrackFinalProperTime(time*s);
	//fastStep.SetPrimaryTrackFinalKineticEnergy(0); //end with no energy
	fastStep.SetPrimaryTrackPathLength(pathLength*mm); //travel calculated distance
	fastStep.SetPrimaryTrackFinalPosition(position); //set to final calculated position
	fastStep.SetTotalEnergyDeposited(ekin_keV*keV); //deposit all energy
	ProcessHit(fastTrack,position,ekin_keV);
	fastStep.KillPrimaryTrack();
	//if(radius <=50){//if it hits the anode radius
	//	G4cout<<"Anode Hit!"<<G4endl;
	//	fSensitive->Hit(fastStep);

	//} else { 
	//	G4cout<<"Electron misses anode"<<G4endl;
	//}
	//ProcessHit(position,ekin_keV,time);
	
	//end
	//fastStep.UpdateStepForPostStep();
	//fastStep.DumpInfo(); //print all step information
	//fastStep.SetKineticEnergy(0);
}

/***** Process secondary electrons created along drift track *****/
void A2HeedModel::ProcessSecondaries(){
	//do nothing for now
	//later deal with secondary electrons created
}


/***** Call a hit in the anode for each electron that reaches it *****/
void A2HeedModel::ProcessHit(const G4FastTrack& fastTrack,G4ThreeVector position, G4double ekin_keV){
//find the volume we are currently working in and set up a touchable
	if (!fNaviSetup)
    {
      fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
      fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
      fNaviSetup = true;
    } else {
      fpNavigator->
        LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
     }
  //fill G4Step with information necessary for the sensitive detector
  //set track to step
  fFakeStep->SetTrack(const_cast<G4Track*>(fastTrack.GetPrimaryTrack()));
  //set touchable for step position
  fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
  //set total energy deposit
  fFakeStep->SetTotalEnergyDeposit(ekin_keV);
	G4VPhysicalVolume* fCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VSensitiveDetector* fSensitive;
	if( fCurrentVolume != 0 ) {
	fSensitive = fCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
	if( fSensitive != 0 ) {
		fSensitive->Hit(fFakeStep);
		//G4cout<<"Hitting sensitive detector"<<G4endl;
		}
	}
}

/***** Create proper response in detector *****/
void A2HeedModel::GenerateDetectorResponse(){
	//do nothing for now
	//later make sure electrons are picked up by anode
}

/***** Reimplement so that class works *****/
void A2HeedModel::ProcessEvent(){
	//reimplement from G4VFastSimulation
}
void A2HeedModel::Reset(){
	//reimplement from G4VFastSimulation
}

