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

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
	G4cout<<"### Detectorconstruction constructor\n";
	checkOverlaps = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials()
{ 
  	//Get nist material manager
  	G4NistManager* nist = G4NistManager::Instance();

    Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	Alluminio = nist->FindOrBuildMaterial("G4_Al");
	Mylar = nist->FindOrBuildMaterial("G4_MYLAR");
	Tungsteno = nist->FindOrBuildMaterial("G4_W");
    Acciaio = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    
  	//G4double densityCF4IC;
  	//G4double densityISOPPAC;
    //densityCF4IC = 0.00014*g/cm3;
  	//densityISOPPAC = 0.00028*g/cm3;           //calcola
  	G4Element* elC = new G4Element("Carbonio", "C", 6., 12.*g/mole);
    G4Element* elH = new G4Element("Idrogeno", "H", 1., 1.*g/mole);
    G4Element* elF =  new G4Element("Fluoro", "F", 9., 18.9984*g/mole);
    
  	IsobutanoPPAC = new G4Material("Isobutano", densityISOPPAC, 2, kStateGas, temperature, pressureISOPPAC);
  	IsobutanoPPAC->AddElement(elC, 4);
    IsobutanoPPAC->AddElement(elH, 10);

    
    tetrafluorometanoIC = new G4Material("tetrafluorometanoIC", densityCF4IC, 2, kStateGas, temperature, pressureCF4IC);
    tetrafluorometanoIC->AddElement(elC, 1);
    tetrafluorometanoIC->AddElement(elF, 4);
	
	//Print all the materials defined
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	G4cout << "### Defining materials in detector construction" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    G4cout << "### Starting detector construction" << G4endl;
    
	//call defined materials
    DefineMaterials();

	//initialization of geometrical variables

   
    //World    
    G4double worldSize = 10.*CLHEP::meter;
    
    G4Box* worldSolid = new G4Box("World",
                                  worldSize/2.,
                                  worldSize/2.,
                                  worldSize/2.);
    
    worldLogic = new G4LogicalVolume(worldSolid,
                                  	Vacuum,
                                    "World");

	worldLogic->SetVisAttributes(G4VisAttributes::Invisible);
	  
    G4PVPlacement* worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "World",
                                                     0,
                                                     false,
                                                     0,
													 checkOverlaps);

	//Geomatry parameters//////////////////////////////////////////////////////////////////////////											

	G4double mylar_PPAC_thickness = 0.9 *um;
	G4double mylar_DC_thickness = 1.5 *um;
	G4double mylar_IC_thickness = 6 *um;

	G4double PPAC_thickness = 120 * mm;
	G4double DC_thickness = 160 * mm;
	G4double DC_spacer_thickness = 120 * mm;
	G4double spacer_front = 15 * mm;
	G4double spacer_back = 35 * mm;
	std::vector<G4double> IC_sections_thickness = {60*mm, 60*mm, 120*mm, 120*mm, 120*mm, 100*mm, 20*mm};

	G4double chamber_width = 1400*mm;
	G4double chamber_height = 500*mm;

	G4double tolerance = 1*mm;

	G4double chamber_thickness = 5*mm;
	std::vector<G4double> thickness = {mylar_PPAC_thickness,
							PPAC_thickness,
							mylar_DC_thickness,
							mylar_DC_thickness,
							DC_thickness,
							DC_spacer_thickness,
							DC_thickness,
							spacer_front,
							mylar_IC_thickness,
							spacer_back,
							IC_sections_thickness[0],
							IC_sections_thickness[1],
							IC_sections_thickness[2],
							IC_sections_thickness[3],
							IC_sections_thickness[4],
							IC_sections_thickness[5],
							IC_sections_thickness[6]};

	std::vector<G4double> cumulative_thickness;


	layers[0] = new Layer(Mylar, mylar_PPAC_thickness, "Mylar_PPAC", false, 0.5, 0.5, 0.5);
	layers[1] = new Layer(IsobutanoPPAC, PPAC_thickness, "Iso_PPAC", false, 0.5, 0.5, 0.5);
	layers[2] = new Layer(Mylar, mylar_DC_thickness, "Mylar_DC0", false, 0.5, 0.5, 0.5);
	layers[3] = new Layer(Mylar, mylar_DC_thickness, "Mylar_DC1", false, 0.5, 0.5, 0.5);
	layers[4] = new Layer(IsobutanoPPAC, DC_thickness, "Iso_DC0", false, 0.5, 0.5, 0.5);
	layers[5] = new Layer(IsobutanoPPAC, DC_spacer_thickness, "Iso_spacer_DC", false, 0.5, 0.5, 0.5);
	layers[6] = new Layer(IsobutanoPPAC, DC_thickness, "Iso_DC1", false, 0.5, 0.5, 0.5);
	layers[7] = new Layer(IsobutanoPPAC, spacer_front, "spacer_front", false, 0.5, 0.5, 0.5);
	layers[8] = new Layer(Mylar, mylar_IC_thickness, "Mylar_IC", false, 0.5, 0.5, 0.5);
	layers[9] = new Layer(tetrafluorometanoIC, spacer_back, "spacer_back", false, 0.5, 0.5, 0.5);
	layers[10] = new Layer(tetrafluorometanoIC, IC_sections_thickness[0], "IC0", false, 0.5, 0.5, 0.5);
	layers[11] = new Layer(tetrafluorometanoIC, IC_sections_thickness[1], "IC1", false, 0.5, 0.5, 0.5);
	layers[12] = new Layer(tetrafluorometanoIC, IC_sections_thickness[2], "IC2", false, 0.5, 0.5, 0.5);
	layers[13] = new Layer(tetrafluorometanoIC, IC_sections_thickness[3], "IC3", false, 0.5, 0.5, 0.5);
	layers[14] = new Layer(tetrafluorometanoIC, IC_sections_thickness[4], "IC4", false, 0.5, 0.5, 0.5);
	layers[15] = new Layer(tetrafluorometanoIC, IC_sections_thickness[5], "IC5", false, 0.5, 0.5, 0.5);
	layers[16] = new Layer(tetrafluorometanoIC, IC_sections_thickness[6], "IC6", false, 0.5, 0.5, 0.5);

	G4double total_length = 0;
	for (const auto & layer: layers) total_length += layer.second->thickness;

	G4double prev_thickness = 0;
	for (const auto & layer: layers){ 
		layer.second->position = prev_thickness +layer.second->thickness/2. - total_length/2.;
		prev_thickness+=layer.second->thickness;
	}

	chamber = new G4Box("chamber", chamber_width/2.+chamber_thickness, chamber_height/2.+chamber_thickness, total_length/2.);
	chamber = new G4SubtractionSolid("chamber2", chamber, new G4Box("chamberdummy",chamber_width/2, chamber_height/2, total_length ));
	chamberLogic = new G4LogicalVolume(chamber, Acciaio, "ChamberLogical");
	G4VisAttributes* chamberVis = new G4VisAttributes(G4Color(0.3, 0.0, 0.2));
	chamberVis->SetForceSolid(true);
	chamberLogic->SetVisAttributes(chamberVis);

	G4ThreeVector shift(0, 0, total_length);

	for (const auto & layer: layers){
		layer.second->solid = new G4Box(layer.second->name.c_str(),
										chamber_width/2.-tolerance,
										chamber_height/2.-tolerance,
										layer.second->thickness/2.);
		layer.second->logical = new G4LogicalVolume(layer.second->solid, 
													layer.second->material,
													std::string(layer.second->name+"Logical").c_str());
		layer.second->logical->SetVisAttributes(layer.second->visualization);
	}

	chamberPhysical = new G4PVPlacement(0,
										shift,
										chamberLogic,
										"chamberPhysical",
										worldLogic,
										false,
										0,
										checkOverlaps);
	
	for (const auto & layer: layers){
		layer.second->physical = new G4PVPlacement(0,
													shift+G4ThreeVector(0, 0, layer.second->position),
													layer.second->logical,
													std::string(layer.second->name+"physical").c_str(),
													worldLogic,
													false,
													0,
													checkOverlaps);
	}

    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	//Sensitive Detector
	G4cout<<"### Setting sensitive detector\n";

	for(const auto & detector: layers){
		if(detector.second->name.compare(0, 2,"IC")==0){
			G4VSensitiveDetector* vDetector = new SensitiveDetector(detector.second->name);
			G4SDManager::GetSDMpointer()->AddNewDetector(vDetector);
			detector.second->logical->SetSensitiveDetector(vDetector);
		}
	}
	G4cout<<"### Finished setting sensitive detector\n";
	
	//Dose Scoring
	//G4MultiFunctionalDetector* multisd = new G4MultiFunctionalDetector("multisd");
	//G4VPrimitiveScorer* dosedet = new G4PSDoseDeposit("dose");
	//multisd->RegisterPrimitive(dosedet);
	//G4SDManager::GetSDMpointer()->AddNewDetector(multisd);
	//SetSensitiveDetector("World",multisd);
		
	//if I use the following expression
	//HPGeLogical->SetSensitiveDetector(multisd);
    //the previois SD does not work anymore (I can associate only one SD to a logical volume)
	//SetSensitiveDetector("HPGeLogical",multisd); 
	//is a workaround!

/*
	G4SDParticleFilter* particleFilter = new G4SDParticleFilter("particleFilter");
	particleFilter->add("e-");  
	dosedet->SetFilter(particleFilter);          
*/   
/*
	G4SDKineticEnergyFilter* energyFilter = new G4SDKineticEnergyFilter("particleFilter");
	energyFilter->SetKineticEnergy(10.*keV,100.*keV); 
	dosedet->SetFilter(energyFilter);          
 */  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..

