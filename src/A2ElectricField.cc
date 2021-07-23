/***** Electric Field Class *****
 * Currently only implements a uniform electric field.
 * Written for use with A2TPC.
 ***** AC Postuma 2021 *****/

#include "A2ElectricField.hh"

#include "G4UniformElectricField.hh" //base class
#include "G4UniformMagField.hh" //magnetic field
#include "G4MagneticField.hh" //also magnetic field - why are you here?
#include "G4FieldManager.hh" //EM field manager
#include "G4TransportationManager.hh" //manages transport in fields
#include "G4EquationOfMotion.hh" //general eqn of motion
#include "G4EqMagElectricField.hh" //eqn of motion in field
#include "G4Mag_UsualEqRhs.hh" //I'm not sure what this does
#include "G4MagIntegratorStepper.hh" //integration of motion in field
#include "G4MagIntegratorDriver.hh" //integration of motion in field
#include "G4ChordFinder.hh" //not sure about this either

#include "G4ExplicitEuler.hh" //different steppers
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/***** Constructor *****/
A2ElectricField::A2ElectricField(){
	//initiate pointers to NULL
	fField=NULL;
	fEquation=NULL;
	fFieldManager=NULL;
	fStepper=NULL;
	fIntegrationDriver=NULL;
	fChordFinder=NULL;
	fStepperType=4; //choose a stepper
	fMinStep=0.010*mm; //initialize minimum step length
}

/***** Destructor ******/
A2ElectricField::~A2ElectricField(){
	//remove things that require manual deletion
	delete fField;
	delete fEquation;
	delete fFieldManager;
	delete fStepper;
	delete fIntegrationDriver;
	delete fChordFinder;
}

/**** Construct field of specified strength *****/
G4ElectricField* A2ElectricField::Construct(G4double fieldStrength){
	fField= new G4UniformElectricField(G4ThreeVector(0.0,0.0,fieldStrength*kilovolt/cm)); //create field
	fEquation= new G4EqMagElectricField(fField); //create equation for field
	fFieldManager = GetGlobalFieldManager(); //get field manager
	UpdateIntegrator(); //set up integrator to do motion in field
	return fField;
}

/***** This function updates integration of motion in field *****/
void A2ElectricField::UpdateIntegrator(){
	CreateStepper(); //pick stepper and get integration driver based on stepper
	fChordFinder = new G4ChordFinder(fIntegrationDriver); //chord finder based on this stepper 
	fFieldManager->SetChordFinder(fChordFinder); //attach to field manager
	fFieldManager->SetDetectorField(fField); //attach to detector
}

/***** This function creates the stepper *****/
void A2ElectricField::CreateStepper(){
	const G4int nvar = 8; //varaibles in each equation
  	auto oldStepper= fStepper;

  	switch ( fStepperType )
  	{
  	  case 0:
  	    fStepper = new G4ExplicitEuler( fEquation, nvar );
  	    G4cout<<"G4ExplicitEuler is calledS"<<G4endl;
  	    break;
  	  case 1:
  	    fStepper = new G4ImplicitEuler( fEquation, nvar );
  	    G4cout<<"G4ImplicitEuler is called"<<G4endl;
  	    break;
  	  case 2:
  	    fStepper = new G4SimpleRunge( fEquation, nvar );
  	    G4cout<<"G4SimpleRunge is called"<<G4endl;
  	    break;
  	  case 3:
  	    fStepper = new G4SimpleHeum( fEquation, nvar );
  	    G4cout<<"G4SimpleHeum is called"<<G4endl;
  	    break;
  	  case 4:
  	    fStepper = new G4ClassicalRK4( fEquation, nvar );
  	    G4cout<<"G4ClassicalRK4 is called"<<G4endl;
  	    break;
  	  case 5:
  	    fStepper = new G4CashKarpRKF45( fEquation, nvar );
  	    G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  	    break;
  	  case 7:
  	    fStepper = 0; // new G4HelixExplicitEuler( fEquation );
  	    G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
  	    break;
  	  case 8:
  	    fStepper = 0; // new G4HelixImplicitEuler( fEquation );
  	    G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
  	    break;
  	  case 9:
  	    fStepper = 0; // new G4HelixSimpleRunge( fEquation );
  	    G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
  	    break;
  	  default:  /* fStepper = 0; // Older code */
  	    fStepper = new G4ClassicalRK4( fEquation, nvar );
  	    G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
  	    break;
  	}

  	delete oldStepper;
  
	//make integration driver based on this stepper
	fIntegrationDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables());
  	//fIntegrationDriver->RenewStepperAndAdjust(fStepper);
}

/****** This function changes the strength of a defined field *****/
void A2ElectricField::SetFieldValue(G4double fieldValue){
	G4ThreeVector fieldVector(0.0,0.0,fieldValue*kilovolt/cm);
	G4FieldManager *fieldManager = GetGlobalFieldManager();
	if (fieldVector != G4ThreeVector(0.,0.,0.)){
		fField = new G4UniformElectricField(fieldVector);
	} else {
		fField=0;
	}
	fieldManager->SetDetectorField(fField);
	fEquation->SetFieldObj(fField);
}

/***** This function is required by the class *****/
void A2ElectricField::GetFieldValue(const G4double point[4], G4double *field) const{
	//reimplement from G4ElectricField
	fField->GetFieldValue(point, field);
}


/***** Also required for this class *****/
G4FieldManager* A2ElectricField::GetGlobalFieldManager(){
	return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
