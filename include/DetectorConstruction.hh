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
// --------------------------------------------------------------
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#endif

#include "G4VUserDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "globals.hh"

#include "map"

#include "SensitiveDetector.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4ExtrudedSolid.hh"
#include "G4Polyhedra.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
//#include "G4PSDoseDeposit.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDParticleFilter.hh"
#include "G4VSDFilter.hh"
#include "G4SDKineticEnergyFilter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:   
    DetectorConstruction();
    ~DetectorConstruction();
    
    void DefineMaterials();

    G4VPhysicalVolume* Construct();
    void ConstructSDandField();

private:
	//varible for checking overlap of volumes
	G4bool checkOverlaps;

    //variabili
    
    //G4double mwppac_pos = 50.*mm;
    G4double temperature = 239; 
    G4double pressureISOPPAC = 0.008*bar; 
    G4double pressureCF4PPAC = 0.008*bar; 
    G4double pressureISOIC = 0.046*bar; 
    G4double pressureCF4IC = 0.046*bar; 
    
	//definition of Materilas
    G4Material* Vacuum;
    G4Material* Alluminio;
    G4Material* Tungsteno;
    G4Material* Acciaio;
	G4Material* Mylar;
	G4Material* LaBr3;
	G4Material* IsobutanoPPAC;
	G4Material* IsobutanoIC;
	G4Material* tetrafluorometanoPPAC;
	G4Material* tetrafluorometanoIC;
	
	//definiton of geometrical variables
    G4int nr_scintillators = 1;

	//definition of Volumes
	G4LogicalVolume* worldLogic;
    
    G4LogicalVolume* LidLogical;
    G4LogicalVolume* ChamberLogical;
    G4LogicalVolume* ChamberDummyLogical;

    std::vector<G4LogicalVolume*> ScintillatorLogical;
    G4LogicalVolume* ScintillatorShellLogical;
    
    //Physical volumes
	G4VPhysicalVolume* LidPhysical;
	G4VPhysicalVolume* ChamberPhysical;
	std::vector<G4VPhysicalVolume*> ScintillatorPhysical;
	G4VPhysicalVolume* ScintillatorShellPhysical;

};


