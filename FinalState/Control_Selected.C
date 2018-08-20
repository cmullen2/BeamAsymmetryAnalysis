{
  //optional Set input file directory, could be done with --hsin or setenv HSIN
  //HSin("/w/work3/home/chris/LatestAnalysisRuns/Data/DataJul18/ChrisOutput/Pi0Analysis/");
  HSin("/w/work3/home/chris/LatestAnalysisRuns/Data/DataJul18/ChrisOutput/Pi0Analysis/");
  //optional set ouput directory or file.
  //Directory=>1 out file for every in file
  //file=> 1 total output file
  //HSout("FULL_PATH_OUTPUT_FILE_DIR_OR_FILE");
  HSout("/w/work3/home/chris/LatestAnalysisRuns/Data/DataJul18/HaspectOutput/Pi0Analysis/AllFilesDev2.0/");
 
  HSfinal("THSProj_Pi0N");

 
  //Set up the input chain to be analysed,could remove HSin() if you want
  TChain* tree=new TChain("HSParticles","datachain");
//  tree->Add(HSin()+"Physics_CBTaggTAPS_14835.root");
  tree->Add(HSin()+"Physics_CBTaggTAPS*");

  //if using proof
  //Give proof destination file (env variables can't be read directly
  if(gProof) {gProof->AddInput(new TNamed("HSOUT",HSout().Data()));gProof->AddInput(new TNamed("HSPROOF","Using Proof"));tree->SetProof();}

  gBenchmark->Start("SelectorTimer");

  tree->Process("Selected.C++");
  
  gBenchmark->Stop("SelectorTimer");
  gBenchmark->Print("SelectorTimer");
  tree->PrintCacheStats();

  //Not actually sure if this is needed or not...
  if(gProof)gProof->ClearData(TProof::kUnregistered | TProof::kForceClear);
}
