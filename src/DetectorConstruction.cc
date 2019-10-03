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
    
  	G4double densityCF4PPAC, densityCF4IC;
  	G4double densityISOPPAC, densityISOIC;
  	densityCF4PPAC = 0.00028*g/cm3;           //calcola
    densityCF4IC = 0.00014*g/cm3;
  	densityISOPPAC = 0.00028*g/cm3;           //calcola
    densityISOIC = 0.00014*g/cm3;
  	G4Element* elC = new G4Element("Carbonio", "C", 6., 12.*g/mole);
    G4Element* elH = new G4Element("Idrogeno", "H", 1., 1.*g/mole);
    G4Element* elF =  new G4Element("Fluoro", "F", 9., 18.9984*g/mole);
    
  	IsobutanoPPAC = new G4Material("Isobutano", densityISOPPAC, 2, kStateGas, temperature, pressureISOPPAC);
  	IsobutanoPPAC->AddElement(elC, 4);
    IsobutanoPPAC->AddElement(elH, 10);

  	IsobutanoIC = new G4Material("Isobutano", densityISOIC, 2, kStateGas, temperature, pressureISOIC);
  	IsobutanoIC->AddElement(elC, 4);
    IsobutanoIC->AddElement(elH, 10);
    
    tetrafluorometanoIC = new G4Material("tetrafluorometanoIC", densityCF4IC, 2, kStateGas, temperature, pressureCF4IC);
    tetrafluorometanoIC->AddElement(elC, 1);
    tetrafluorometanoIC->AddElement(elF, 4);

    tetrafluorometanoPPAC = new G4Material("tetrafluorometanoIC", densityCF4PPAC, 2, kStateGas, temperature, pressureCF4IC);
    tetrafluorometanoPPAC->AddElement(elC, 1);
    tetrafluorometanoPPAC->AddElement(elF, 4);
	
	//LaBr3
 	G4Element* La = new G4Element("Lanthanum", "La", 57.,138.90547*g/mole);
  	G4Element* Br = new G4Element("Bromium", "Br", 35.,79.904*g/mole);	

	LaBr3 = new G4Material("LaBr3",  5.07 * g / cm3,  2);
	LaBr3->AddElement(La, 1);
	LaBr3->AddElement(Br, 3);

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
	G4double LidThickness = 5*mm;
	G4double LidRadius = 407*mm;
	G4double ChamberThickness = 25*mm;
	G4double ChamberHeight = 745*mm;
	G4double LidHeight = 231*mm;
	G4double LidTheta = atan(LidRadius/LidHeight);
	G4double Lid3DRadius =LidRadius/sin(LidTheta); 
	G4double DummyDist =300*mm;
	G4double tmpDist =30*mm;
	//Chamber Body////////////////////////////////////////////////////////////////////////////////////////////////
	G4double Tolerance = 0*mm;
/*  	G4double ChamberR[9]={0*mm, LidRadius, LidRadius, LidRadius+45*mm, LidRadius+45*mm, 
						LidRadius-25*mm, LidRadius-25*mm, 0*mm, 0*mm};
	G4double ChamberZ[9]={-(ChamberHeight-(Lid3DRadius-LidHeight))-ChamberThickness/2.-Tolerance, -(ChamberHeight-(Lid3DRadius-LidHeight))-ChamberThickness/2.-Tolerance,
						Lid3DRadius-LidHeight-25*mm-ChamberThickness/2.-Tolerance, Lid3DRadius-LidHeight-25*mm-ChamberThickness/2.-Tolerance, Lid3DRadius-LidHeight-ChamberThickness/2.-Tolerance,
						Lid3DRadius-LidHeight-ChamberThickness/2.-Tolerance,-(ChamberHeight-(Lid3DRadius-LidHeight))+25*mm-ChamberThickness/2.-Tolerance,
						 -(ChamberHeight-(Lid3DRadius-LidHeight))+25*mm-ChamberThickness/2.-Tolerance,-(ChamberHeight-(Lid3DRadius-LidHeight))-ChamberThickness/2.-Tolerance};
 	G4VSolid* ChamberBody = new G4Polycone("ChamberBodyUnfinished",
										0*deg, 360*deg,
										9, ChamberR, ChamberZ);
 */

	G4VSolid* ChamberBody = new G4Tubs("ChamberBody1",LidRadius-25*mm,LidRadius, ChamberHeight,0*deg, 360*deg);
	
	ChamberBody = new G4IntersectionSolid("ChamberBody2",
										ChamberBody,new G4Box("Boxtmp",
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										ChamberHeight/2.),
										0,G4ThreeVector(-DummyDist, -DummyDist, 0*mm));

	ChamberBody = new G4DisplacedSolid("ChamberBody3", ChamberBody, 0,G4ThreeVector(0,0,-ChamberHeight/2.+Lid3DRadius-LidHeight-ChamberThickness/2.));
	ChamberBody = new G4UnionSolid("ChamberBody4", ChamberBody, new G4Tubs("TubsTmp", 0, LidRadius,ChamberThickness, 0*deg, 360*deg), 
									0, G4ThreeVector(0,0,-(ChamberHeight-(Lid3DRadius-LidHeight))-ChamberThickness-ChamberThickness/2.));

	ChamberBody = new G4IntersectionSolid("ChamberBody5",
										ChamberBody,new G4Box("Boxtmp",
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										ChamberHeight),
										0,G4ThreeVector(-DummyDist, -DummyDist, 0*mm));

	G4VSolid* ChamberDummy = new G4Box("ChamberDummy", 
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										ChamberHeight);

	ChamberDummy = new G4DisplacedSolid("ChamberDummy2", ChamberDummy, 0, G4ThreeVector(-DummyDist+ChamberThickness, -DummyDist+ChamberThickness, 0*mm));
	ChamberDummy = new G4SubtractionSolid("ChamberDummy3", ChamberDummy, 
										new G4Box("ChamberDummy4", 
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										ChamberHeight+ChamberThickness), 0 ,G4ThreeVector(-DummyDist, -DummyDist, 0*mm));
	ChamberDummy = new G4IntersectionSolid("ChamberDummy5", ChamberDummy, 
											new G4Tubs("ChamberDummy6", 0, LidRadius, ChamberHeight/2+ChamberThickness, 0*deg, 360*deg), 
											0, 
											G4ThreeVector(0,0,-ChamberHeight/2.+Lid3DRadius-LidHeight-ChamberThickness-0.5*ChamberThickness));
/*										
	ChamberDummy = new G4SubtractionSolid("ChamberDummy3",
										new G4Box("Boxtmp",
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										ChamberHeight),
										ChamberDummy,
										0,G4ThreeVector(+DummyDist-ChamberThickness, +DummyDist-ChamberThickness, 0*mm));
*/


/*	G4VSolid* ChamberDummy2 = new G4Box("ChamnberDummy2", 
										DummyDist+188*mm+190*mm-tmpD,
										DummyDist+170*mm+166*mm-tmpDist,
										Lid3DRadius);
	ChamberDummy = new G4SubtractionSolid("ChamberDummy3", ChamberDummy2, ChamberDummy);
	ChamberBody = new G4UnionSolid("ChamberBody6", ChamberBody, ChamberDummy);
*/

	//Cahmber Lid///////////////////////////////////////////////////////////////////////////////////////

	G4VSolid* ChamberLidDummy1 = new G4Sphere("ChamberLidDummy1",
										Lid3DRadius-LidThickness,Lid3DRadius, 
										0.*deg, 360*deg,
										0.*deg, LidTheta);
	G4VSolid* ChamberLidDummy2 = new G4Box("ChamberLidDummy2",
										DummyDist+188*mm+190*mm-tmpDist,
										DummyDist+170*mm+166*mm-tmpDist,
										Lid3DRadius);
	G4VSolid* ChamberLidDummy3 = new G4Tubs("ChamberLidDummy3",
								LidRadius-LidThickness,LidRadius+45*mm, 25*mm/2., 0*deg, 200*deg);
	G4VSolid* ChamberLid = new G4IntersectionSolid("ChamberLid",ChamberLidDummy1, ChamberLidDummy2, 
												0,G4ThreeVector(-DummyDist, -DummyDist, 0*mm));

	//Valve
	G4double ValveR[9]={50*mm, 50*mm+LidThickness, 50*mm+LidThickness, 50*mm+LidThickness+20*mm, 
						50*mm+LidThickness+20*mm, 0*mm, 0*mm, 50*mm, 50*mm};
	G4double ValveZ[9]={0*mm, 0*mm, Lid3DRadius, Lid3DRadius, Lid3DRadius+50*mm, Lid3DRadius+50*mm,
						Lid3DRadius+17.5*mm, Lid3DRadius+17.5*mm, 0};

	G4VSolid* ValveDummy = new G4Polycone("ValveDummy",
										0*deg, 360*deg,
										9, ValveR, ValveZ);
	
	ValveDummy = new G4DisplacedSolid("ValveDummyDisplaced", ValveDummy, 0, G4ThreeVector(0*mm, 200*mm, 0*mm));

	G4VSolid* Valve =new G4SubtractionSolid("Valve", ValveDummy, new G4Sphere("SphereDummy",
										0,Lid3DRadius-LidThickness, 
										0.*deg, 360*deg,
										0.*deg, 180*deg));
							
	//Obj
	DummyDist = 300*mm;
	G4VSolid* ObjDummy = new G4Box("ObjDummy", LidRadius-DummyDist/2., ChamberThickness/2., Lid3DRadius-32*mm-DummyDist/2.);
	ObjDummy = new G4DisplacedSolid("ObjDummyDisplaced", ObjDummy, 0,G4ThreeVector(-DummyDist/2., 50*mm,DummyDist/2.));
	G4RotationMatrix* ObjRotation = new G4RotationMatrix();
	ObjRotation->rotateZ(-30*deg);
	ObjDummy = new G4DisplacedSolid("ObjDummyDisplaced2", ObjDummy, ObjRotation,G4ThreeVector(0*mm, 0*mm, 0*mm));


	G4VSolid* ObjDummy2 = new G4Box("ObjDummy2", LidRadius-DummyDist/2., ChamberThickness/2., Lid3DRadius-32*mm-DummyDist/2.);
	ObjDummy2 = new G4DisplacedSolid("ObjDummyDisplaced", ObjDummy2, 0,G4ThreeVector(-DummyDist/2., -50*mm,DummyDist/2.));
	ObjDummy2 = new G4DisplacedSolid("ObjDummyDisplaced2", ObjDummy2, ObjRotation,G4ThreeVector(0*mm, 0*mm, 0*mm));

	ObjDummy = new G4UnionSolid("ObjDummyUnited", ObjDummy, ObjDummy2);

	G4VSolid* ObjDummy3 = new G4Box("ObjDummy3",160*mm/2.,175*mm/2., ChamberThickness/2.);
	ObjDummy3 = new G4DisplacedSolid("ObjDummy3Displaced", ObjDummy3, 0, G4ThreeVector(0*mm,LidRadius-175*mm/2,Lid3DRadius-32*mm+ChamberThickness/2.));
	ObjRotation->rotateZ(-90*deg);
	ObjDummy3 = new G4DisplacedSolid("ObjDummy3Displaced2", ObjDummy3, ObjRotation, G4ThreeVector(0, 0, 0));

	G4VSolid* Obj = new G4SubtractionSolid("Obj", ObjDummy,new G4Sphere("SphereDummy",
										0,Lid3DRadius-LidThickness, 
										0.*deg, 360*deg,
										0.*deg, 180*deg));
							
	Obj = new G4UnionSolid("Obj2", Obj, ObjDummy3);
	Obj = new G4SubtractionSolid("Obj2", Obj, ChamberBody);
	//Flange
	G4RotationMatrix* RotationFlange = new G4RotationMatrix();
	RotationFlange->rotateZ(230 * deg);

	//Unions =)
	ChamberLid = new G4UnionSolid("ChamberLidNoValve", ChamberLid, ChamberLidDummy3, RotationFlange,G4ThreeVector(0,0,Lid3DRadius-LidHeight));

	ChamberLid = new G4UnionSolid("ChamberLidNoObj", ChamberLid, Valve);

	ChamberLid = new G4UnionSolid("ChamberLid", ChamberLid, Obj);


	//Scintillator/////////////////////////////////////////////////////////////////////////////////////////

	G4double ScintillatorRadius = 38.1*mm;
	G4double ScintillatorLength = 76.2*mm;
	G4double ShellThichkness = 5*mm;
	G4double PMTLength = 177*mm;
	G4double ToleranceScintillator = 2*mm;
	G4double ShellRadius = 82*mm/2;


	G4VSolid* Scintillator = new G4Tubs("Scintillator", 0, ScintillatorRadius,ScintillatorLength/2., 0*deg, 360*deg);
	
	//ScintillatorShell////////////////////////////////////////////////////////////////////////////////


	G4VSolid* ScintillatorShell = new G4Tubs("ScintillatorShell1",ScintillatorRadius, 
												ShellRadius, (ScintillatorLength+PMTLength)/2., 0*deg, 360*deg );
	G4VSolid* ScintillatorFront = new G4Tubs("ScintillatorFront", 0, ShellRadius, ShellThichkness/2., 0*deg, 360*deg);
	ScintillatorFront = new G4DisplacedSolid("ScintillatorFront", ScintillatorFront, 0, G4ThreeVector(0, 0, -(ScintillatorLength+ShellThichkness)/2));
	ScintillatorShell = new G4DisplacedSolid("ScintillatorShell", ScintillatorShell, 0,G4ThreeVector(0, 0, PMTLength/2));
	ScintillatorShell = new G4UnionSolid("ScintillatorShell", ScintillatorShell, ScintillatorFront);

	//Logical////////////////////////////////////////////////////////////////////////////////////////////////	
	LidLogical = new G4LogicalVolume(ChamberLid, Acciaio, "ChamberLidLogical");
	G4VisAttributes* ChamberLidAttribute = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
	ChamberLidAttribute->SetForceSolid(true);
	LidLogical->SetVisAttributes(ChamberLidAttribute);


	ChamberLogical = new G4LogicalVolume(ChamberBody, Acciaio, "ChamberBodyLogical");
	G4VisAttributes* ChamberAttribute = new G4VisAttributes(G4Color(0.3, 0.3, 0.3));
	ChamberAttribute->SetForceSolid(true);
	ChamberLogical->SetVisAttributes(ChamberAttribute);

	ChamberDummyLogical = new G4LogicalVolume(ChamberDummy, Acciaio, "ChamberDummyLogical");
	G4VisAttributes* ChamberDummyAttribute = new G4VisAttributes(G4Color(0.3, 0.3, 0.3));
	ChamberDummyAttribute->SetForceSolid(true);
	ChamberDummyLogical->SetVisAttributes(ChamberDummyAttribute);


	ScintillatorShellLogical = new G4LogicalVolume(ScintillatorShell, Alluminio, "ScintillatorShellLogical");
	G4VisAttributes* ScintillatorShellAttribute = new G4VisAttributes(G4Color(0., 1, 0.));
	ScintillatorShellAttribute->SetForceSolid(true);
	ScintillatorShellLogical->SetVisAttributes(ScintillatorShellAttribute);

	//Physical////////////////////////////////////////////////////////////////////////////////////////////////////

	G4RotationMatrix* LidRotation = new G4RotationMatrix();
	//LidRotation->rotateX(90*deg);
	//LidRotation->rotateZ(180*deg);
	LidPhysical = new G4PVPlacement(LidRotation, 
									G4ThreeVector(0, 0, 0),
									LidLogical,
									"LidPhysical",
									worldLogic,
									false,
									0,
									checkOverlaps);
 
	LidPhysical = new G4PVPlacement(LidRotation, 
									G4ThreeVector(0, 0, 0),
									ChamberLogical,
									"ChamberPhysical",
									worldLogic,
									false,
									0,
									checkOverlaps);

	new G4PVPlacement(LidRotation, 
									G4ThreeVector(0, 0, 0),
									ChamberDummyLogical,
									"ChamberDummyPhysical",
									worldLogic,
									false,
									0,
									checkOverlaps);
  


	G4double ScintillatorDistLid =18*mm;
	G4double ScintillatorDistEach = 2*ScintillatorRadius+2*ShellThichkness+3*cm;//The choice is on 10*cm
	G4double ScintillatorDistEach2 = 2*ScintillatorRadius+2*ShellThichkness+2*cm;//The choice is on 10*cm
	G4double ThetaMax = LidTheta*0.8;
	G4int NrRings = ThetaMax*(Lid3DRadius+ScintillatorLength)/ScintillatorDistEach2;




	for (int jj=0; jj<NrRings; jj++){
	G4double ThetaRing = ThetaMax-jj*(ThetaMax/NrRings);
	G4int nrRing = 2*pi*sin(ThetaRing)*(Lid3DRadius+ScintillatorLength)/ScintillatorDistEach;

	std::map <int,std::map<int, bool>> Disable;
	for (int pp=0; pp<nrRing; pp++){
		Disable[pp][0]=true;
	}
	Disable[0][1]=true;
	Disable[0][2]=true;
//	Disable[0][3]=true;
	Disable[1][1]=true;
	Disable[1][2]=true;
//	Disable[1][3]=true;


	Disable[17][1]=true;
	Disable[11][2]=true;
//	Disable[5][3]=true;

	Disable[5][1]=true;
	Disable[6][1]=true;
	Disable[7][1]=true;


	Disable[3][2]=true;
	Disable[4][2]=true;
	Disable[5][2]=true;
	//da togliere
/*	Disable[7][1]=true;
	Disable[5][1]=true;
	Disable[12][0]=true;
	Disable[8][0]=true;
	Disable[2][3]=true;
*/	

	for (int ii=0; ii<nrRing;ii++ ){
			G4cout<<"LidRadius :"<<Lid3DRadius/mm<<std::endl;
			if (Disable[ii][jj]) continue;
			ScintillatorLogical.push_back( new G4LogicalVolume(Scintillator, LaBr3, "ScintillatorLogical"));
			G4VisAttributes* ScintillatorAttribute = new G4VisAttributes(G4Color(1, 0., 0.));
			ScintillatorAttribute->SetForceSolid(true);
			ScintillatorLogical.back()->SetVisAttributes(ScintillatorAttribute);
			// 	ScintillatorShellPhysical = new G4PVPlacement(Ring1Rotation,G4ThreeVector(0,0,Lid3DRadius+ScintillatorDistLid).rotateX(-ThetaMax),
			G4ThreeVector tmpVec(0, 0, Lid3DRadius + ScintillatorDistLid + ChamberThickness + LidThickness);
			tmpVec.rotateX(-ThetaRing);
			G4double offset=0;
			if (jj==3) offset= 0.5* 360 *deg /nrRing;
			tmpVec.rotateZ(ii * 360 * deg / nrRing+offset);

			G4RotationMatrix* tmpRot = new G4RotationMatrix();
			tmpRot->rotate(ThetaRing,tmpVec.cross(G4ThreeVector(0,0,1)));	
						new G4PVPlacement(tmpRot, tmpVec,
								 ScintillatorShellLogical,
								 "ScintillatorShellPhysical",
								 worldLogic,
								 false,
								 0,
								 checkOverlaps);
 
			G4VPhysicalVolume* tmpScintillatorPhysical = new G4PVPlacement(tmpRot, tmpVec,
									ScintillatorLogical.back(),
									"ScintillatorPhysical",
									worldLogic,
									false,
									0,
									checkOverlaps);
			ScintillatorPhysical.push_back(tmpScintillatorPhysical);
								
	}
	}
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	//Sensitive Detector
	G4cout<<"### Setting sensitive detector\n";
	G4int ii=0;
	for(const auto & Detector: ScintillatorLogical){
		std::string det_name("det");
		det_name.append(std::to_string(ii));
//		det_name.append(std::to_string(ii));
//		G4cout << det_name <<"\n";
		G4VSensitiveDetector* vDetector = new SensitiveDetector(det_name);
		G4SDManager::GetSDMpointer()->AddNewDetector(vDetector);
		Detector->SetSensitiveDetector(vDetector);
		ii++;
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

