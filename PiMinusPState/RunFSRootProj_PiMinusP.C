//root --hsdata --hsfinal=THSProj_PiMinusP RunFSRootProj_PiMinusP.C
//You need to replace Proj_PiMinusP with your final state class name
{
  //Create FinalState
  THSProj_PiMinusP* fs=new THSProj_PiMinusP();
  // fs->SetGenerated(); //just analyse generated branch
   fs->SetMaxParticles(100);
  //create datamanager
  THSDataManager* dm=new THSDataManager();
  TChain chain("HSParticles");
  chain.Add("/w/work1/home/chris/GenVectMethod/ChargedPions/Phys*.root");
  dm->InitChain(&chain);
  //connect Project to HSParticles
  fs->SetDataManager(dm);
  Int_t counter=0;
  
  //create ouput tree
  TFile* outfile=new TFile("/w/work1/home/chris/GenVectMethod/ChargedPions/HaspectOut/Dev1.0Files20Physics.root","recreate");
  TTree* outtree=new TTree("FinalTree","output tree");
  fs->FinalStateOutTree(outtree); //connect ouput tree to project branches
  
  gBenchmark->Start("timer");
  
  while(dm->ReadEvent()){//loop over events
    fs->ProcessEvent();
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  outfile->cd();
  outtree->Write();
  delete outfile;
}
