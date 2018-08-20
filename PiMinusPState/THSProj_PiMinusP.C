/**
 * \class THSProj_PiMinusP
 * 
 * Template final class .
 * 
 * Users should create their own analysis specific project classes.
 * 
 */



#include "TDatabasePDG.h"
#include "THSProj_PiMinusP.h"
#include <algorithm>


THSProj_PiMinusP::THSProj_PiMinusP(){
  SetVerbose(1);
  // CheckCombi(); //comment out to remove messages
  
  //Set final state detected particles
  //AddParticle(particle,true/false you want to write in final vector, genID for linking to generated truth value)
  //Note if particle is added to final with a valid genID it will be used
  //to determine the correct permutation of the simulated event
  AddParticle(&fPhoton,kTRUE,3);
  AddParticle(&fPim,kTRUE,1);
  AddParticle(&fProton,kTRUE,2);
  AddParticle(&fProtonSpec,kTRUE,0);

  //Set final state parents
  
  
  
  TString PID("ALL"); //set this to which particles you want to id via pdg code alone, or set it in individual AddTopology
  TString INCLUSIVE("Beam");//set this to which particles you want reaction to be inclusive of, or set it in individual AddTopology "ALL"=> completely inclusive

  //include topology for analysis and get index
  AddTopology("Beam:pi-:proton",
             bind(&THSProj_PiMinusP::Init_Iter0, this),
             bind(&THSProj_PiMinusP::Topo_0, this),
             PID,INCLUSIVE);

  
  
  THSFinalState::InitFS();
}

void THSProj_PiMinusP::FileStart(){
  //May be called when a new file is opened
  
}

//Define topology Iterator functions
void THSProj_PiMinusP::Init_Iter0(){
  //THSParticle iterator initialisation
  //For topology pi-:proton

   AutoIter(); //Let finalstate try and work out the iterattor for you, you can remove this if you want to do it yourself
}

//Define topology functions
void THSProj_PiMinusP::Topo_0(){
  //For topology pi-:proton
  fCorrectedProtonEnergy = ProtonELossCorrection(fProton.P4p()->Theta(), fProton.P4p()->E());
  fMissNucleon = fPhoton.P4() + fFreeProton - fPim.P4();

}

void THSProj_PiMinusP::Kinematics(){
  //Do calculations if Good Event

 if(fProton.P4p()->Phi()==0 && fProton.P4p()->Theta()==0){
    fDetErrs = -1;
  }
  else{
    fDetErrs = 0;
  }

  //Timing
  fTagTime = fPhoton.Time();
  fPimTime = fPim.Time();

  //Adjustment for jitter, common start and stop different in Taps vs Ball
  if( fPim.Detector()<8)  fPimTagDiffTime = fPimTime - fTagTime;
  if( fPim.Detector()>7)  fPimTagDiffTime = -fPimTime - fTagTime;




//Apply Mikhail's energy correction to the Proton so can use in Doms function to calculate pim 4vec
  fProton.P4p()->SetE(fCorrectedProtonEnergy + fProton.PDGMass() );

//Decision to be made on which energy to take, calc proton from pim or calc pim from proton
  fProtonCalcEnergy = CalcQFThreeBodyRecoilPartT(fPhoton.P4().E(), fPim.P4(), fProton.P4(),massTarget,massProton,massProton); 


//Check Dom's function worked correctly
  if(fProtonCalcEnergy<0){
    fDomFuncErrs = -1;
  }
  else{
    fDomFuncErrs = 0;
  }

  if(std::isnan(fProtonCalcEnergy)){
    fEnergyErrs = -1;

  }
  else{
    fEnergyErrs = 0;
    fProton.P4p()->SetE(fProtonCalcEnergy + fProton.PDGMass() );

  }


  if(fDetErrs==0) fProton.TakePDGMassFromE();

/*  fPimCalcEnergy = CalcQFThreeBodyRecoilPartT(fPhoton.P4.E(), fProton.P4(), fPim.P4(),massTarget,massPim,massProton); //Proton spectator for this channel
  if(fPimCalcEnergy<0){
    fDomFuncErrs = -1;
  }
  else{
    fDomFuncErrs = 0;
  }

  if(std::isnan(fPimCalcEnergy)){
    fEnergyErrs = -1; 
    
  }
  else{
    fEnergyErrs = 0;
    fPim.P4p()->SetE(fPimCalcEnergy + fPim.PDGMass() );

  }
  
  if(fDetErrs==0) fProton.TakePDGMassFromE();
  fPim.TakePDGMassFromE();

*/

//  fNucleonEnergyFinal = fProton.P4p()->E();
//  fPimEnergyFinal = fPim.P4p()->E();
  
  //Reconstruct missing or combined particles
  fBeamEnergy = fPhoton.P4p()->E();
  fSpectator = fPhoton.P4() + fTarget - fProton.P4() - fPim.P4();
  fPimMassDiff = (TDatabasePDG::Instance()->GetParticle(-211)->Mass()) *1000 - fPim.P4p()->M();
  fMissP4 = fPhoton.P4() + fTarget - fProton.P4();  // now using corrected proton fPim.P4();
  fMissMass = fMissP4.M();
  fInvMass = fPim.P4p()->M();
  fSpecMass = fSpectator.M(); //Should be almost junk
  fSpecMom = fSpectator.P();
  fCoplanarity = ( (ROOT::Math::VectorUtil::DeltaPhi(fPim.P4(), -(fProton.P4()) ) )*TMath::RadToDeg());
  fInitialN = fTarget - fSpectator;
  fConeAngle =ROOT::Math::VectorUtil::Angle(fMissNucleon,fProton.P4().Vect());
 
	  
  //Detector Tests
  fProtonDetector = fProton.Detector();
  fPimDetector = fPim.Detector();

  //Polarisation
  fTaggChannel = fPhoton.Detector();
  fPolState = fEventInfo->TarPolDir();
  fPol = fEventInfo->TarPol();
  fPolStateD = fEventInfo->TarPolDir();
 // fPolStateD = 1; // For SIMS


  //Zero Polarisation
  if(fPol==0){
    fPolErrs = -1;
  }
  else{
    fPolErrs=0;
  }

 //Using Kinematics functions
  fKine.SetGammaTarget(fPhoton.P4(),fInitialN);
  fKine.SetMesonBaryon(fPim.P4(),fProton.P4());
  fKine.PhotoCMDecay();
  fCMPhi = (fKine.Phi())*TMath::RadToDeg();
  fCosth = fKine.CosTheta(); 
  fW = fKine.W();
  fProduc = fFreeProton + fPhoton.P4();
  fWII = fW - fProduc.M();

 if(fDetErrs==-1){     //No longer need energyerrs or dom   fDetErrs==-1 || fEnergyErrs==-1 || fDomFuncErrs==-1 ){
    fAnyErrs = -1;
  }
  else{
    fAnyErrs=0;
  }

}
//////////////////////////////////////////////////
/// Define conditions for track to be considered
/// good in event. Adding conditions on junk tracks
///  can greatly reduce combitorials etc.
/// kFALSE=> track ignored completely
Bool_t THSProj_PiMinusP::CheckParticle(THSParticle* part){
  return kTRUE;
}

void THSProj_PiMinusP::FinalStateOutTree(TTree* tree){
  THSFinalState::fFinalTree=tree;
  //tree->Branch("Final",&fFinal);//If you want to save the final THSParticles
  //  tree->Branch("MissMass",&fMissMass,"MissMass/D");
  tree->Branch("Topo",&fTopoID,"Topo/I");
  tree->Branch("NPerm",&fNPerm,"NPerm/I");
  tree->Branch("NDet",&fNDet,"NDet/I");

  tree->Branch("MissMass",&fMissMass,"MissMass/D");
  tree->Branch("Phi",&fCMPhi,"Phi/D");   //COM Phi angle of meson
  tree->Branch("Costh",&fCosth,"Costh/D");  //COM costheta angle of meson 
  tree->Branch("SpecMom",&fSpecMom,"SpecMom/D"); //Spectator mass
  tree->Branch("Coplanarity",&fCoplanarity,"Coplanarity/D"); //Coplanarity between pi0 and part. nucl.
  tree->Branch("BeamEnergy",&fBeamEnergy,"BeamEnergy/D");
  tree->Branch("W",&fW,"W/D");
  tree->Branch("InvMass",&fInvMass,"InvMass/D"); //Inv mass of meson
  tree->Branch("Pol",&fPol,"Pol/D"); //Polarisation mag. also has state encompassed
  tree->Branch("TaggChannel",&fTaggChannel,"TaggChannel/D"); //Tagger Channels
  tree->Branch("ConeAngle",&fConeAngle,"ConeAngle/D");
  tree->Branch("DetErrs",&fDetErrs,"DetErrs/D"); // Error indicator for NO CB or TAPS hit
  tree->Branch("PimTagDiffTime",&fPimTagDiffTime,"PimTagDiffTime/D");
  tree->Branch("AnyErrs",&fAnyErrs,"AnyErrs/D"); // Any Error given from the above three 
  tree->Branch("Correct",&fCorrect,"Correct/I");
  tree->Branch("DCorrect",&fDCorrect,"DCorrect/D"); //Does GenSim match TopoSim
  tree->Branch("SpecMass",&fSpecMass,"SpecMass/D"); //Spectator mass
  tree->Branch("WII",&fWII,"WII/D");
  tree->Branch("CorrectedProtonEnergy",&fCorrectedProtonEnergy,"CorrectedProtonEnergy/D");
  tree->Branch("TagTime",&fTagTime,"TagTime/D");
  tree->Branch("PimTime",&fPimTime,"PimTime/D");
  //  tree->Branch("NucleonEnergyFinal",&fNucleonEnergyFinal,"NucleonEnergyFinal/D");
  tree->Branch("PimMassDiff",&fPimMassDiff,"PimMassDiff/D");
  tree->Branch("ProtonDetector",&fProtonDetector,"ProtonDetector/D");
  tree->Branch("PimDetector",&fPimDetector,"PimDetector/D");
  tree->Branch("PolState",&fPolState,"PolState/I");
  tree->Branch("PolStateD",&fPolStateD,"PolStateD/D");
  tree->Branch("PolErrs",&fPolErrs,"PolErrs/D");





}
