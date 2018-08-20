{
  THSSkeleton* sk=new THSSkeleton();
  //Give your project a name
  sk->SetFinalState("Proj_PiMinusP",kTRUE);  //creating the project class with perms
  sk->SetFinalStateTopo("pi-:proton");
  //Set the actual final state particles of the reaction
  sk->SetFinalStateParts("Photon:gamma,Pim:pi-,Proton:proton,ProtonSpec:proton");
  sk->CreateMyFinalState(); 
}
