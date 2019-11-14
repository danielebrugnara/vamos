//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "EventAction.hh"
#include "SensitiveDetectorHit.hh"
#include "Analysis.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction(RunAction* runAction):
fSensitiveDetector_ID(-1),
fVerboseLevel(10),
fRunAction(runAction),
Edep(0.)
{
	G4cout << "### EventAction instantiated" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*)
{
	//G4cout << "### Begin of event" << G4endl;
	Edep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* aEvent)
{   
	//instantiating The Sensitive Detector Manager
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

    //instantiating The Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

//	std::vector<G4double> hitsE;
//	std::vector<G4int> hitsDet;
        
	for (int ii=0; ii<7; ii++){
	    //Hit Detection System
		if(fSensitiveDetector_ID == -1||1) {
	    	G4String SensitiveDetectorName="IC";
			SensitiveDetectorName.append(std::to_string(ii));
	        if(SDman->FindSensitiveDetector(SensitiveDetectorName,0)) {
				SensitiveDetectorName.append("/collection");
	            fSensitiveDetector_ID = SDman->GetCollectionID(SensitiveDetectorName);
//				G4cout << "Able to find detector "<<SensitiveDetectorName<<std::endl;
	        }else{
				continue;
//				G4cout << "Unable to find detector "<<SensitiveDetectorName<<std::endl;
			}
	    }else{
	//		G4cout <<"Sensitive det ID = "<< fSensitiveDetector_ID<<std::endl;
		}

	    SensitiveDetectorHitsCollection* fSensitiveDetectorHC = 0;
	    G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();

	    if(HCE) {
	        if(fSensitiveDetector_ID != -1) {
	            G4VHitsCollection* aHC = HCE->GetHC(fSensitiveDetector_ID);
	            fSensitiveDetectorHC = (SensitiveDetectorHitsCollection*)(aHC);
	        }
	    }


	//	G4double ti = 0.;
	//	G4double tf = 0.;
	//	G4double tHit = 0.;

		//filling the scoring ntuple
		G4double eHit = 0.;
		G4double Edet=0;
	    if(fSensitiveDetectorHC ) {
	        int vNumberOfHit = fSensitiveDetectorHC->entries();
	        for(int i=0; i<vNumberOfHit; i++) {
	            SensitiveDetectorHit* aHit = (*fSensitiveDetectorHC)[i];
//				tHit = aHit->GetTime();
				eHit = aHit->GetEnergy();
				if(i==0){
					G4double tmpX = aHit->GetPos().getX();
					G4double tmpY = aHit->GetPos().getY();
					G4double tmpZ = aHit->GetPos().getZ();

					analysisManager->FillNtupleDColumn(0,7+0,tmpX/CLHEP::mm);
					analysisManager->FillNtupleDColumn(0,7+1,tmpY/CLHEP::mm);
					analysisManager->FillNtupleDColumn(0,7+2,tmpZ/CLHEP::mm);
				}
//				Edep += eHit;
				Edet += eHit;
//				if (i==0) {
//					ti = tHit;  
//					tf = tHit;  
//				}
//				if (tHit < ti && eHit > 0.) ti = tHit;
//				if (tHit > tf && eHit > 0.) tf = tHit;
	        }
//			G4cout <<"i value :" <<ii<< "    Edet: "<<Edet<<"\n";
//			hitsE.push_back(Edep/CLHEP::keV);
//			hitsDet.push_back(ii);
		}
	//	G4cout << "Edep " << Edep << std::endl;
		analysisManager->FillNtupleDColumn(0,ii,Edet/CLHEP::MeV);
	}
//	if (1) {
//		for (int ii=0; ii<50; ii++){
//		   	analysisManager->FillNtupleDColumn(0,0,hitsE);
//		}	
//	   	analysisManager->FillNtupleDColumn(0,1,ti/CLHEP::ns);
//	   	analysisManager->FillNtupleDColumn(0,2,tf/CLHEP::ns);
		//analysisManager->FillNtupleDColumn(0,3,eHit/CLHEP::keV);
	analysisManager->AddNtupleRow(0);   
//	}
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

