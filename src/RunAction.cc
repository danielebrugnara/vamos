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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"
#include "Analysis.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4SDManager.hh"

#include "G4AccumulableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():
G4UserRunAction(),
fMessenger(0)
{
	G4RunManager::GetRunManager()->SetPrintProgress(1000);
    
    fMessenger = new RunActionMessenger(this);

    fFileName = "SimuOutput";
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "### Analysis Manager is " << analysisManager->GetType() << " ###" << G4endl;
    
    //Creating the scoring ntuple
    analysisManager->CreateNtuple("ICDet","ICDet");
	for (int ii=0; ii<7; ii++){ 
		std::string name = "IC";
		name.append(std::to_string(ii));
    	analysisManager->CreateNtupleDColumn(name);
	}
    analysisManager->CreateNtupleDColumn("PosX");
    analysisManager->CreateNtupleDColumn("PosY");
    analysisManager->CreateNtupleDColumn("PosZ");
    analysisManager->CreateNtupleDColumn("E0");
//    analysisManager->CreateNtupleDColumn("ti");
//    analysisManager->CreateNtupleDColumn("tf");
    analysisManager->FinishNtuple();
      	
  	G4cout << "### RunAction instantiated" << G4endl;    	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun() {return new Run;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  	G4cout << "### Begin of RunAction" << G4endl;    	
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fFileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{  	   
    //write ntuples (results from SD and SteppingAction) on a file
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();

    
    G4int nofEvents = run->GetNumberOfEvent();
  	if (nofEvents == 0) return;
  
  	//results from the Scorer
  	const Run* aRun = static_cast<const Run*>(run);
  	G4int nbGoodEvents = aRun->GetNbGoodEvents(); 
  	//G4double sumDose   = aRun->GetSumDose();                           
        
  	//print
  	if (IsMaster()) {
    	G4cout
     	<< G4endl
     	<< "--------------------End of Global Run-----------------------"
     	<< G4endl
     	<< " The run had " << nofEvents << " events";
  	} else {
    	G4cout
     	<< G4endl
     	<< "--------------------End of Local Run------------------------"
     	<< G4endl
     	<< " The run had " << nofEvents << " events";
  	}      
  	
  	G4cout
    << "; Nb of 'good' events: " << nbGoodEvents << G4endl 
    << G4endl
    //<< " Total dose in solid: " << sumDose/gray << " Gray" << G4endl  
    << "------------------------------------------------------------" << G4endl 
    << G4endl;


	//Write Results from Scorer on a text file  
    if (IsMaster()) {
    	std::ofstream fFileOut;
        fFileOut.open("SimulationSummary.dat", std::ofstream::out | std::ofstream::app);
        fFileOut << "seed: " << G4Random::getTheSeed() 
        		// << " | " << sumDose/gray << " Gray in HPGe"    
        	 	 << std::endl;
        fFileOut.close();
  	}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetFileName(G4String filename){if(filename != "") {fFileName = filename;}}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

