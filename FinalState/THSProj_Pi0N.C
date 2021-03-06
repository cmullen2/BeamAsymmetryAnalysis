/**
 * \class THSProj_Pi0N
 * 
 * Template final class .
 * 
 * Users should create their own analysis specific project classes.
 * 
 */



#include "TDatabasePDG.h"
#include "THSProj_Pi0N.h"
#include <algorithm>
#include "THSTopology.h"

THSProj_Pi0N::THSProj_Pi0N(){
  SetVerbose(1);
 // CheckCombi(); //comment out to remove messages
  
  //Set final state detected particles
  //AddParticle(particle,true/false you want to write in final vector, genID for linking to generated truth value)
  //Note if particle is added to final with a valid genID it will be used
  //to determine the correct permutation of the simulated event
  AddParticle(&fPhoton,kTRUE,4);
  AddParticle(&fProton,kTRUE,1);
  AddParticle(&fGamma3,kTRUE,1); //Call this gamma3 and keep neutron in topo function make sure has right mass etc. and link to 3rd gamma
  AddParticle(&fGamma1,kTRUE,2);
  AddParticle(&fGamma2,kTRUE,3);

  //Set final state parents
  AddParticle(&fPion,kTRUE,-1);
  ConfigParent(&fPion,&fGamma1);
  ConfigParent(&fPion,&fGamma2);
  
  
  
  TString PID("ALL"); //set this to which particles you want to id via pdg code alone, or set it in individual AddTopology (Proton or None can test,BEAM:Proton)
  TString INCLUSIVE("Beam");//set this to which particles you want reaction to be inclusive of, or set it in individual AddTopology "ALL"=> completely inclusive (Beam)

  //include topology for analysis and get index
  AddTopology("Beam:gamma:gamma:gamma",
  	      bind(&THSProj_Pi0N::Init_Iter0, this),
  	      bind(&THSProj_Pi0N::Topo_0, this),
  	      PID,INCLUSIVE);

   AddTopology("Beam:proton:gamma:gamma",
   	      bind(&THSProj_Pi0N::Init_Iter1, this), 
   	      bind(&THSProj_Pi0N::Topo_1, this),
   	      PID,INCLUSIVE);
   //The proton topology contains the third gamma, need to watch out for it using it in the actual code, Should not iterate it and should keep it under gamma3 
  
  
  THSFinalState::InitFS();
}

void THSProj_Pi0N::FileStart(){
  //May be called when a new file is opened
  
}

//Define topology Iterator functions
void THSProj_Pi0N::Init_Iter0(){
  //THSParticle iterator initialisation
  //For topology gamma:gamma:gamma


  //THSParticle iterator initialisation

  THSParticleIter* diter=CreateParticleIter((fMapPDGtoParticle[-22]),1);
  AddSelectToSelected(diter,1,1,&fPhoton);
  diter=CreateParticleIter(&fVecGams,3);
  THSParticleIter* pi0diter=AddSelectToSelected(diter,1,2,&fGamma1,&fGamma2);
  AddSelectToRemainder(pi0diter,1,1,&fGamma3);

  //If create iterator must configure it
  // fDetIter[fTID_1]->ConfigureIters();
  //  AutoIter(); //Let finalstate try and work out the iterattor for you, you can remove this if you want to do it yourself

}
void THSProj_Pi0N::Init_Iter1(){
  //THSParticle iterator initialisation
  //For topology proton:gamma:gamma

  //THSParticle iterator initialisation
  THSParticleIter* diter=CreateParticleIter((fMapPDGtoParticle[-22]),1);
  AddSelectToSelected(diter,1,1,&fPhoton);
  diter=CreateParticleIter(&fVecGams,2);
  AddSelectToSelected(diter,1,2,&fGamma1,&fGamma2);
  diter=CreateParticleIter(&fVecProtons,1);
  AddSelectToSelected(diter,1,1,&fProton);

  //If create iterator must configure it
  //  fDetIter[fTID_0]->ConfigureIters();
  //  AutoIter(); //Let finalstate try and work out the iterattor for you, you can remove this if you want to do it yourself
}

//Define topology functions
void THSProj_Pi0N::Topo_0(){
  //For topology gamma:gamma:gamma

  fTopologies=-1; //Neutron Channel
  fPion.SetP4(fGamma1.P4() + fGamma2.P4());
  fPion.SetTime((fGamma1.Time() + fGamma2.Time() )/2 );

  //A dummy representative of Mikhail's energy correction.(kinetic energy of the neutron)
  fCorrectedProtonLabEnergyb4 =fNeutron.P4p()->E(); //If used to setE then have this as neutron.GetE etc.
  
  fMissNucleon=fPhoton.P4() + fFreeNeutron - fPion.P4();

  //Dominik's function
  fNucleonLabEnergy = CalcQFThreeBodyRecoilPartT(fPhoton.P4().E(), fPion.P4(), fNeutron.P4(), massTarget, massNeutron,massProton);
 
  //Set up a generic nucleon variable to be used in kinematics
  //fNeutron.SetPDGcode(2112);
  fNeutron.SetP4(fGamma3.P4());
  fNeutron.SetTime(fGamma3.Time());
  fNucleon = &fNeutron;

}
void THSProj_Pi0N::Topo_1(){
  //For topology proton:gamma:gamma

  fTopologies=1; //Proton Channel
  fPion.SetP4(fGamma1.P4() + fGamma2.P4());
  fPion.SetTime((fGamma1.Time() + fGamma2.Time() )/2 );
  // Mikhails energy correction function
  fCorrectedProtonLabEnergyb4 = ProtonELossCorrection(fProton.P4p()->Theta(), fProton.P4p()->E());
  fMissNucleon=fPhoton.P4() + fFreeProton - fPion.P4();
  //Dominik's function
  fNucleonLabEnergy = CalcQFThreeBodyRecoilPartT(fPhoton.P4().E(), fPion.P4(), fProton.P4(), massTarget, massProton,massNeutron);
  //Set up a generic nucleon variable to be used in kinematics
  fNucleon = &fProton;

}


void THSProj_Pi0N::Kinematics(){
  //Do calculations if Good Event


  //  Test for ball or taps hit Make this a detector based selection
  if(fNucleon->P4p()->Phi()==0 && fNucleon->P4p()->Theta()==0){
    fDetErrs= -1;
  }
  else{
    fDetErrs=0;
  }
  
  // Test for error in dominiks energy function
  if(fNucleonLabEnergy<0){
    fDomFuncErrs=-1;
  }
  else{
    fDomFuncErrs=0;
  }

  // Timing
  fTagTime=fPhoton.Time();
  //Adjustment for jitter, common start/stop
  if(  fGamma1.Detector()<8 && fGamma2.Detector()<8)fGammaAveTagDiffTime= fPion.Time() - fTagTime;
  if(  fGamma1.Detector()>7 && fGamma2.Detector()>7)fGammaAveTagDiffTime= -fPion.Time() - fTagTime;
  if(  fGamma1.Detector()>7 && fGamma2.Detector()<8)fGammaAveTagDiffTime= -fGamma1.Time() - fTagTime;
  if(  fGamma1.Detector()<8 && fGamma2.Detector()>7)fGammaAveTagDiffTime= -fGamma2.Time() - fTagTime;


  if(std::isnan(fNucleonLabEnergy) ){
    fEnergyErrs =  -1;
    fNucleon->P4p()->SetE(fNucleon->P4p()->E()+fNucleon->PDGMass());  //Not sure if 1000 needed see if still in GeV.
  }
  else{
    fEnergyErrs = 0;
    fNucleon->P4p()->SetE(fNucleonLabEnergy + fNucleon->PDGMass() ); 
  }

  if(fDetErrs==0) fNucleon->TakePDGMassFromE(); //Hits with no cb or taps will cause vector stretch error if not dealt with by fDetErrs

  //Reconstruct missing or combined particles
  fBeamEnergy=fPhoton.P4p()->E();
  fSpectator=fPhoton.P4() +  fTarget - fNucleon->P4() -fPion.P4(); //Is this correct????
  fM2gamma= (TDatabasePDG::Instance()->GetParticle(111)->Mass())*1000  - fPion.P4p()->M(); //in MeV now //MassDiff(); //Need to convert to MeV the pdg mass
  fMissP4=fPhoton.P4() + fTarget - fPion.P4();
  fMissMass=fMissP4.M();
  fMissPionP4 = fPhoton.P4() +fTarget -fNucleon->P4();
  fMissMassPion = fMissPionP4.M();
  fInvMass=fPion.P4p()->M();   //MeasMass();
  fSpecMass=fSpectator.M();
  fSpecMom=fSpectator.P();
  fCoplanarity=(ROOT::Math::VectorUtil::DeltaPhi( fPion.P4(), -(fNucleon->P4()) )  )*TMath::RadToDeg();
  fInitialN=fTarget- fSpectator;
  fConeAngle=ROOT::Math::VectorUtil::Angle(fMissNucleon,fNucleon->P4().Vect());
 


  // //Polarisation
   fTaggChannel = fPhoton.Detector();
   fPolState=fEventInfo->TarPolDir();
   fPol=fEventInfo->TarPol();   //Edge plane =1,-1,0  describes Para,Perp or Moeller( not necess in that order)
   //Double for Fitting in Roofit as a binvar
   fPolStateD =fEventInfo->TarPolDir(); //For Prod Data
   fPolStateD = -1; // For  Sim Data so that have both polarisation for binning for fitting  (Need to rerun for +1) DO the same for polstate
   fPolState = -1; //SIMS


   //Zero Polarisation 
   if (fPol==0){
     fPolErrs=-1; 
   }
   else{
     fPolErrs=0;
   }

  //Detectors used for each particle
  fDetector = fNucleon->Detector();
  // Centre of mass Calculations
  fKine.SetGammaTarget(fPhoton.P4(),fInitialN);
  fKine.SetMesonBaryon(fPion.P4(),fNucleon->P4());
  fKine.PhotoCMDecay();
  fCMPhi=(fKine.Phi())*TMath::RadToDeg() ;
  fCosth=fKine.CosTheta(); 
  fW = fKine.W();
  fProduc=fFreeProton + fPhoton.P4() ;
  fWII =fW - fProduc.M() ;

  if(fDetErrs==-1 || fEnergyErrs==-1 || fDomFuncErrs==-1  ){
    fAnyErrs=-1;
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
Bool_t THSProj_Pi0N::CheckParticle(THSParticle* part){
  return kTRUE;
}

void THSProj_Pi0N::FinalStateOutTree(TTree* tree){
  THSFinalState::fFinalTree=tree;
  //tree->Branch("Final",&fFinal);//If you want to save the final THSParticles
  //  tree->Branch("MissMass",&fMissMass,"MissMass/D");
  tree->Branch("Topo",&fTopoID,"Topo/I");
  tree->Branch("NPerm",&fNPerm,"NPerm/I");
  tree->Branch("NDet",&fNDet,"NDet/I");

  tree->Branch("TopoMine",&fTopologies,"TopoMine/D");//Two topologies proton particapant and neutron part,pi0 
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
  tree->Branch("EnergyErrs",&fEnergyErrs,"EnergyErrs/D");  // Error indicator for Dom's func.
  tree->Branch("DetErrs",&fDetErrs,"DetErrs/D"); // Error indicator for NO CB or TAPS hit
  tree->Branch("GammaAveTagDiffTime",&fGammaAveTagDiffTime,"GammaAveTagDiffTime/D"); 
  tree->Branch("DomFuncErrs",&fDomFuncErrs,"DomFuncErrs/D"); // Negative Energy from Doms Function indicator 
  tree->Branch("AnyErrs",&fAnyErrs,"AnyErrs/D"); // Any Error given from the above three 
  tree->Branch("MissMass",&fMissMass,"MissMass/D"); //Missing Mass of Participant
  tree->Branch("MissMassPion",&fMissMassPion,"MissMassPion/D"); //Missing Mass of Pion
  tree->Branch("Correct",&fCorrect,"Correct/I"); 
    tree->Branch("DCorrect",&fDCorrect,"DCorrect/D"); //Does GenSim match TopoSim
  tree->Branch("SpecMass",&fSpecMass,"SpecMass/D"); //Spectator mass
  tree->Branch("WII",&fWII,"WII/D"); 
  
  tree->Branch("PolState",&fPolState,"PolState/I"); //PolState +-1 or 0 ,Perp Para moeller
  tree->Branch("PolStateD",&fPolStateD,"PolStateD/D"); //PolState +-1 or 0 ,Perp Para moeller
  tree->Branch("PolErrs",&fPolErrs,"PolErrs/D"); // Any Error given from the above three 
 
}



/*
Bool_t THSProj_Pi0N::WorkOnEvent(){
  //Should this event be saved?
  THSFinalState::fGoodEvent=kTRUE;
  THSFinalState:: fCorrect=0; //Correct permutation? used for simulation only
  //If generated MC events
  if(THSFinalState::fIsGenerated) Init_Generated();
  else{//Look for reconstructed events
    if(FindInclusiveTopology(-22)==-1){fGoodEvent=kFALSE;return fIsPermutating0=kFALSE;} 

    //Do they correspond to a defined topology?
    else if(fCurrTopo==fTID_0) Topo_0();
    else if(fCurrTopo==fTID_1) Topo_1();
    else fGoodEvent=kFALSE;
    //Get truth values
    Init_Generated();
  }


  //Calc kinematics
  Kinematics();
if(fGammaAveTagDiffTime<-100 || fGammaAveTagDiffTime>100){
fGoodEvent=kFALSE;

}
  
  //Check if assigned vectors agree with true generated
  //Simulation only
  THSFinalState::CheckTruth();
  //SIMS NEXT TWO LINES NEEDED BUT CRASH DATA FILES
  fGamma1.SetTruth(frGenParts->at(3));
  fGamma2.SetTruth(frGenParts->at(2));  //See fConfig class with parent.

  THSFinalState::CheckTruth();


  //For Sims Temporary fix!
  fDCorrect=fCorrect -1;
  //For Prod, Emp (Fix this at some stage! Needed so that they give the same value for cuts during fitting)
//  fDCorrect=fCorrect;

  if(fIsGenerated) return kTRUE; //Generated only 1 permutation
  return kTRUE;
}
*/


Bool_t THSProj_Pi0N::WorkOnEvent(){
  //Should this event be saved?
  THSFinalState::fGoodEvent=kTRUE;
  THSFinalState::fCorrect=0; //Correct permutation? used for simulation only
  //If generated MC events
  InitGenerated();
  if(!THSFinalState::fIsGenerated){
    //Look for reconstructed events
    //if reconstructed Find the detected particles in this event
    //Found a topology execute its Topo function
      fCurrTopo->Exec();
  }
  
  //Calc kinematics
  Kinematics();
  if(fGammaAveTagDiffTime<-100 || fGammaAveTagDiffTime>100){
    fGoodEvent=kFALSE;
  }
  
  //Check if assigned vectors agree with true generated
  //Simulation only
  THSFinalState::CheckTruth();
  fGamma1.SetTruth(frGenParts->at(3));
  fGamma2.SetTruth(frGenParts->at(2));  //See fConfig class with parent.
  THSFinalState::CheckTruth();

  //For Sims Temporary fix!
  fDCorrect=fCorrect -1;
  //For Prod, Emp (Needed so that they give the same value for cuts during fitting)
//  fDCorrect=fCorrect;

  if(THSFinalState::fIsGenerated) return kTRUE; //Generated only 1 permutation
  return kTRUE;
}

