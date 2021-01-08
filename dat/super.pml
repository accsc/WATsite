load prot_amber.pdb
load MyProtein.pdb
super prot_amber, MyProtein
remove solvent
remove resname Na\+ or resname Cl-
save ref.pdb, prot_amber
quit

