
#ifndef THSPROJ_PIMINUSP_h
#define THSPROJ_PIMINUSP_h

#include "THSFinalState.h"
#include "THSParticle.h"
#include <vector>

class THSProj_PiMinusP : public THSFinalState{

 public :
  THSProj_PiMinusP();
  ~THSProj_PiMinusP()=default;


  //Init functions
  void Init_Iter0();
   //Topology action functions
  void Topo_0();
  void Init_Generated();
  //void Init_IterX();
  //void Topo_X();
  virtual void FileStart();
  virtual Bool_t  CheckParticle(THSParticle* part);

  void Kinematics();
  protected :

  
  //Initial state
  HSLorentzVector fTarget=HSLorentzVector(0,0,0,1875.612);

  //Final Particles Detected
  THSParticle fPhoton=THSParticle(-22);
  THSParticle fPim=THSParticle("pi-");
  THSParticle fProton=THSParticle("proton");
  THSParticle fProtonSpec=THSParticle("proton");

  //For Kinematics
  HSLorentzVector fMissNucleon;
  HSLorentzVector fFreeProton=HSLorentzVector(0,0,0,938.272);
  HSLorentzVector fSpectator;
  HSLorentzVector fMissP4;
  HSLorentzVector fProduc;
  HSLorentzVector fInitialN;



  //Discriminators
  Double_t fCorrectedProtonEnergy=0;
  Double_t fDetErrs=0;
  Double_t fDomFuncErrs=0;
  Double_t fTagTime=0;
  Double_t fPimTime=0;
  Double_t fPimTagDiffTime=0;
  Double_t fEnergyErrs=0;
  Double_t fProtonCalcEnergy=0;
  Double_t fPimCalcEnergy=0;
  Double_t fBeamEnergy=0;
  Double_t fPimMassDiff=0;
  Double_t fMissMass=0;
  Double_t fInvMass=0;
  Double_t fSpecMass=0;
  Double_t fSpecMom=0;
  Double_t fCoplanarity=0;
  Double_t fConeAngle=0;
  Double_t fConePhi=0;
  Double_t fProtonDetector=0;
  Double_t fPimDetector=0;
  Double_t fTaggChannel=0;
  Double_t fPol=0;
  Double_t fPolStateD=0;
  Double_t fPolErrs=0;
  Double_t fCMPhi=0;
  Double_t fW=0;
  Double_t fCosth=0;
  Double_t fWII=0;
  Double_t fAnyErrs=0;
  Double_t fDCorrect=0;
  Int_t fPolState=0;


  //For QFEnergy Function
  Double_t massTarget=1875.612;
  Double_t massProton=938.272; 
  Double_t massPim=139.57;


   public :
  virtual void FinalStateOutTree(TTree* tree);

  //______________________________________________________________________________
  Double_t ProtonELossCorrection(Double_t ProtTheta, Double_t NaIDeposit)
  {
    //Mikhails proton energy correction for the polarimeter Mainz beamtime of Aug 16
    // Returns kinetic energy of proton from the deposit in the sodium iodide and proton theta angle(in radians?)  (test case of theta=90deg or pi/2 rads gives no change in Ekin)

    Double_t E_kin=NaIDeposit+(201.915-57.9314*sin(ProtTheta))*(exp((-0.000800067-0.00451967*sin(ProtTheta))*NaIDeposit))+(-82.3023+23.2409*sin(ProtTheta));

    return (E_kin);

  }


  //______________________________________________________________________________
  Double_t CalcQFThreeBodyRecoilPartT(Double_t beamE, HSLorentzVector p4Meson, HSLorentzVector p4Part,
				      Double_t mTarg, Double_t mPart, Double_t mSpec)
  {
    // Calculate the kinetic energy of the recoil participant with mass 'mPart'
    // in meson photoproduction off a target with mass 'mTarg' with beam energy
    // 'beamE' in a quasi-free three-body-decay:
    //
    // gamma + Target(Participant + Spectator) -> Participant + Spectator + Meson
    //
    // Use the reconstructed theta and phi angles in the 4-vector of the
    // participant 'p4Part' and the fully reconstructed 4-vector
    // of the meson 'p4Meson' to calculate the kinetic energy of the recoil
    // participant. The spectator mass is 'mSpec'.

    // set input kinematics variables
    Double_t mesonPx = p4Meson.Px();
    Double_t mesonPy = p4Meson.Py();
    Double_t mesonPz = p4Meson.Pz();
    Double_t mesonE  = p4Meson.E();

    Double_t partTheta = p4Part.Theta();
    Double_t partPhi   = p4Part.Phi();

    // calculate terms
    Double_t a = mesonE - beamE - mTarg;
    Double_t b = (mesonE + mPart - beamE - mTarg) * (mesonE + mPart - beamE - mTarg)
      - mSpec * mSpec
      - mesonPx * mesonPx
      - mesonPy * mesonPy
      - mesonPz * mesonPz
      - beamE * beamE
      + 2 * beamE * mesonPz;
    Double_t c =   mesonPx * TMath::Sin(partTheta) * TMath::Cos(partPhi)
      + mesonPy * TMath::Sin(partTheta) * TMath::Sin(partPhi)
      + (mesonPz - beamE ) * TMath::Cos(partTheta);

    Double_t quadEquA = a*a - c*c;
    Double_t quadEquB = a*b - 2*c*c*mPart;
    Double_t quadEquC = b*b / 4;

    return (-quadEquB + TMath::Sqrt(quadEquB*quadEquB - 4*quadEquA*quadEquC)) / (2*quadEquA);
  }

 
};

#endif //ifdef THSProj_PiMinusP
