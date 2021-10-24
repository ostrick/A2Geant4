//ACP 2021
//implementation of Heed processes in Geant4 to code electron drift in TPC
//based off method described in https://arxiv.org/pdf/1806.05880.pdf
//and implemented in https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE

#ifndef A2DriftModel_h
#define A2DriftModel_h 1

#include "G4VFastSimulationModel.hh"

class A2Target;
class A2SD;

class A2DriftModel : public G4VFastSimulationModel {
public:
	A2DriftModel(G4String, G4Region*, A2Target*,A2SD*); //constructor
	~A2DriftModel(); //destructor

	virtual G4bool IsApplicable(const G4ParticleDefinition&); //return true for electron
	virtual G4bool ModelTrigger(const G4FastTrack&); //trigger if particle energy is below 5 keV
	virtual void DoIt(const G4FastTrack&, G4FastStep&); //pass control of particle to the model

	virtual void ProcessEvent(); //record relevant data
	virtual void Reset(); //reset class variables

protected:
	virtual void Transport(G4FastStep&, const G4FastTrack&, G4String, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double); //move the electron through the active volume
	void ProcessHit(const G4FastTrack&,G4ThreeVector,G4double,G4double);
	void SetConstants(G4Region*);
	G4double drift_vel; //drift velocity
	G4double long_diff; //longitudinal diffusion
	G4double trans_diff; //transverse diffusion
	G4double fHePressure; //from TPC file -  to pick correct data
	G4int fHeIsotope;
	A2Target* fA2Target;
	A2SD* fA2SD;
	G4Step* fFakeStep;
	G4TouchableHandle fTouchableHandle;
	G4Navigator* fpNavigator;
	G4bool fNaviSetup;
	G4StepPoint* fFakePreStepPoint;
        G4StepPoint* fFakePostStepPoint;
	G4Track* fFakeTrack;
};

#endif
