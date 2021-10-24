/***** Electron drift model for TPC *****
 * Implementation of TPC physics simulated in Heed/Garfield/Degrad.
 * Based off method described in https://arxiv.org/pdf/1806.05880.pdf and
 * https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE.
 * Uses G4VFastSimulation to implement manual definition of electron drift.
 * Gaussians based on work by Fabian Metzger, presentation shared by PNPI.
 ***** AC Postuma 2021 *****/

#include "A2DriftModel.hh"

#include "G4VPhysicalVolume.hh" //volumes
#include "G4Electron.hh" //particle the model applies to
#include "G4SystemOfUnits.hh" //units
#include "A2DetectorConstruction.hh" //detector construction
#include "G4RunManager.hh" //run
#include "A2Target.hh" //targets
#include "A2SD.hh" //sensitive detector
#include "G4TrackingManager.hh" //tracking of particles
#include "G4EventManager.hh" //events
#include "G4TransportationManager.hh" //particle transport
#include "G4VVisManager.hh" //track visualization
#include "CLHEP/Random/RandGauss.h" //random generation

using namespace CLHEP;

/**** Constructor *****/
A2DriftModel::A2DriftModel(G4String modelName, G4Region* actVol, A2Target* target, A2SD* anode)
: G4VFastSimulationModel(modelName, actVol), fA2Target(target), fA2SD(anode) { //fast simulation implements user defined physics response
	//initiate pointers
	fFakeStep = new G4Step(); //step used to call hit in SD
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint(); //step point
  	fFakePostStepPoint = fFakeStep->GetPostStepPoint(); //step point
  	fTouchableHandle   = new G4TouchableHistory(); //touchable for step
  	fpNavigator        = new G4Navigator(); //navigator to find SD
	fNaviSetup = false; //if setup has already been done
	//read in data to set drift velocity, diffusion coefficients
	SetConstants(actVol);
}

/***** Destructor *****/
A2DriftModel::~A2DriftModel(){
	//remove objects requiring manual deletion
	delete fFakeStep;
  	delete fpNavigator;
}

/***** Called in SteppingAction: checks particle type if model is applicable ******/
G4bool A2DriftModel::IsApplicable(const G4ParticleDefinition& particleType){
	G4String particleName = particleType.GetParticleName();
	if(particleName=="e-")return true; //only applicable to electrons
	return false;
}

/***** Called in SteppingAction: conditions in which to trigger model *****/
G4bool A2DriftModel::ModelTrigger(const G4FastTrack& fastTrack){
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	if (ekin <=1)return true; //trigger for kinetic energy below 1 keV
	return false;
}

/***** This function contains the main operation of the model *****/
void A2DriftModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep){
	//get all relevant step data from the track
	G4ThreeVector direction = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
	G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition()/mm;
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
	G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	//make sure the electron is in a position to be transported
	G4bool position = true; //it normally will be
	if (worldPosition.z() < -115.5) position = false; //electron already past anode
	G4double radius = sqrt(worldPosition.x()*worldPosition.x()+worldPosition.y()*worldPosition.y());
	if ((worldPosition.z() > 115.5)&&(radius<5)) position=false; //electron behind cathode
	if (position == true){
		Transport(fastStep, fastTrack, particleName, ekin, time, worldPosition.x(), worldPosition.y(), worldPosition.z(), direction.x(), direction.y(), direction.z());
	}
	}

//Generate electron position, time when reaching anode
void A2DriftModel::Transport(G4FastStep& fastStep,const G4FastTrack& fastTrack, G4String particleName, double ekin_keV, double t, double x_mm, double y_mm, double z_mm, double dx, double dy, double dz){
	//G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	/****transport each electron to anode ****/
	G4double z_pos = -115.5; //set final z position to anode z position
	G4double pathLength = z_pos - z_mm; //total length to final z position in mm

	//use Gaussians as defined by Fabian Metzger in his TPC work
	//set means and sigmas of Gaussian distribtuions
	G4double mean_x=x_mm; //mean for Gaussian calc of x pos
	G4double mean_y=y_mm; //mean for Gaussian calc of y pos
	G4double mean_t = abs(pathLength/drift_vel); //mean for Gaussian calc of time: comes out in ms
	G4double sigma_diff = trans_diff*sqrt(abs(pathLength)); //sigma for Gaussian calc of x,y positions: close enough to mm
	G4double sigma_time = long_diff/drift_vel*sqrt(abs(pathLength)); //comes out close enough to ms

	//use random number generation to get values for positions and times
	RanluxEngine aRandEngine; //make an engine for random number generation
	RandGauss gaussian(aRandEngine);
	G4double x_pos = gaussian.shoot(mean_x,sigma_diff); //calculate an x position: mm
	G4double y_pos = gaussian.shoot(mean_y,sigma_diff); //calc y mm
	G4double time = gaussian.shoot(mean_t,sigma_time); //calc a time ms

	//combine position data into a vector
	G4ThreeVector position = G4ThreeVector(x_pos*mm,y_pos*mm,z_pos*mm);

	//let's see what's up here
	//G4cout<<"("<<x_mm<<","<<y_mm<<","<<z_mm<<") to ("<<x_pos<<","<<y_pos<<","<<z_pos<<")"<<G4endl;
	//G4cout<<"r="<<sqrt(x_mm*x_mm+y_mm*y_mm)<<"mm to r="<<sqrt(x_pos*x_pos+y_pos*y_pos)<<"mm in "<<time<<" ms"<<G4endl;
	
	/**** set final track parameters *****/
	fastStep.SetPrimaryTrackFinalProperTime(time);
	fastStep.SetPrimaryTrackPathLength(pathLength*mm); //travel calculated distance
	fastStep.SetPrimaryTrackFinalPosition(position); //final calculated position
	fastStep.SetTotalEnergyDeposited(ekin_keV*keV); //deposit all energy
	
	/**** kill step and call hit ****/
	ProcessHit(fastTrack,position,ekin_keV,time);
	fastStep.KillPrimaryTrack();
}


/***** Call a hit in the anode for each electron that reaches it *****/
void A2DriftModel::ProcessHit(const G4FastTrack& fastTrack,G4ThreeVector position, G4double ekin_keV, G4double drift_time){
/**** set up touchable in current volume ****/
	if (!fNaviSetup) {
		fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
		fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
		fNaviSetup = true;
    	} else {
      		fpNavigator->
        	LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
     	}

  	/**** fill G4Step with information necessary for the sensitive detector ****/
  	//set track to step
  	fFakeStep->SetTrack(const_cast<G4Track*>(fastTrack.GetPrimaryTrack()));
	//set touchable for step position
  	fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
	fFakePreStepPoint->SetPosition(position); //try to get it to work properly
  	//G4cout<<position<<" "<<sqrt(position.x()*position.x()+position.y()*position.y())<<G4endl;
	//set total energy deposit
  	fFakeStep->SetTotalEnergyDeposit(ekin_keV);
  	//set the time of hit
  	fFakeStep->GetPreStepPoint()->SetGlobalTime(drift_time*ms); 

	/**** call hit in sensitive detector ****/
  	G4VPhysicalVolume* fCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VSensitiveDetector* fSensitive;
	if( fCurrentVolume != 0 ) { //if you can find volume
	fSensitive = fCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
	if( fSensitive != 0 ) { //if volume has sensitive detector
		fSensitive->Hit(fFakeStep); //call hit for the fake step
		//G4cout<<"Calling new hit"<<G4endl;
		}
	}
}

/***** Reimplement so that class works *****/
void A2DriftModel::ProcessEvent(){
	//reimplement from G4VFastSimulation
}
void A2DriftModel::Reset(){
	//reimplement from G4VFastSimulation
}

/**** Assign gas parameters depending on isotope, pressure of helium ****/
void A2DriftModel::SetConstants(G4Region *gasRegion){
	G4String name=gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetName();
	G4double pressure = gasRegion->GetRootLogicalVolumeIterator()[0]->GetMaterial()->GetPressure()/bar;
	if (name.contains("3")){
		drift_vel=2840+61*(30-pressure); //linear approx: real data from 30 bar
		trans_diff=0.0325+0.00072*(30-pressure);
		long_diff=0.0214+0.00049*(30-pressure);
	} else {
		drift_vel=2630+55*(30-pressure); //linear approx: real data from 30 bar
		trans_diff=0.0353+0.00059*(30-pressure);
		long_diff=0.0219+0.00071*(30-pressure);
	}
	G4cout<<name<<" "<<pressure<<" "<<drift_vel<<" "<<trans_diff<<" "<<long_diff<<G4endl;
}
