//root --hsdata --hsfinal=THSProj_Pi0N RunFSLundProj_Pi0N.C
//You need to replace Proj_Pi0N with your final state class name
{
  //Create FinalState
  THSProj_Pi0N* fs=new THSProj_Pi0N();
  //create datamanager
  THSLundReader* dm=new THSLundReader();
  dm->Init("/lund/file/name.lund","HSParticles");
  fs->SetDataManager(dm);
  Int_t counter=0;
  
  //create ouput tree
  TFile* outfile=new TFile("gentest.root","recreate");
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
