{
  THSSkeleton* sk=new THSSkeleton();
  //Give your project a name
  sk->SetFinalState("Proj_Pi0N",kTRUE);  //creating the project class with perms
  //Set the detected particle combinations you will analyse
//  sk->SetProjectTopo("pi+:pi-,pi+:pi-:pi0,proton:pi+:pi-:pi0");
  sk->SetFinalStateTopo("gamma:gamma:gamma,proton:gamma:gamma");
  //Set the actual final state particles of the reaction
  //These will just be used as the data member names
 // sk->SetProjectFinal("Electron:e-,Proton:proton,Pip:pi+,Pim:pi-");
//  sk->SetFinalStateParts("Photon:gamma,Proton:proton,Neutron:neutron,Pion:pi0,Gamma1:gamma,Gamma2:gamma");
  sk->SetFinalStateParts("Photon:gamma,Proton:proton,Neutron:neutron,Gamma1:gamma,Gamma2:gamma");
  sk->SetFinalStateParents("Pion:pi0;Gamma1;Gamma2");
  //Make some code
  sk->CreateMyFinalState(); 
}
