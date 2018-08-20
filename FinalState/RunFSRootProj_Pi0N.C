//root --hsdata --hsfinal=THSProj_Pi0N RunFSRootProj_Pi0N.C
//You need to replace Proj_Pi0N with your final state class name
{
  //Create FinalState
  THSProj_Pi0N* fs=new THSProj_Pi0N();
  // fs->SetGenerated(); //just analyse generated branch
   fs->SetMaxParticles(250);
  //create datamanager
  THSDataManager* dm=new THSDataManager();
  TChain chain("HSParticles");
//  chain.Add("/indir/out_*root");
  chain.Add("/w/work1/home/chris/GenVectMethod/NeutralPions/TestNew/Physics_CBTaggTAPS_14946.root");
//  chain.Add("/w/work3/home/chris/LatestAnalysisRuns/Data/DataJul18/ChrisOutput/Pi0Analysis/Physics_CBTaggTAPS_1*");
  dm->InitChain(&chain);
  //connect Project to HSParticles
  fs->SetDataManager(dm);
  Int_t counter=0;
  
  //create ouput tree
  TFile* outfile=new TFile("/w/work1/home/chris/GenVectMethod/NeutralPions/HaspectOut/TestDev2.0Files1PhysicsV2.root","recreate");
//  TFile* outfile=new TFile("/w/work3/home/chris/LatestAnalysisRuns/Data/DataJul18/HaspectOut/Pi0Analysis/Dev1.0FilesAllPhysics.root","recreate");
  TTree* outtree=new TTree("HSParticles","output tree");
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
