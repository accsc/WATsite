#!/bin/tcsh

setenv WATSITEHOME   /programs/WATsite/
setenv mybindingsite 'MyBindingSite.pdb' #this is the ligand or residue that deendifnes the binding site
setenv myprotein     'MyProtein.pdb'

# choose from (0, 1, 2]
# 0 = run reduce
# 1 = rename reisdues (HIE/HID/HIP, LYS/LYN, ASP/ASH, GLU/GLH) based on pre-generated hydrogen atoms by e.g. propka...
# 2 = amberFF compatible residue names 
setenv protonation   0 

# choose from (0, 1]
# 0 = no ligand in the binding site
# 1 = ligand is simulated together in the protein 
setenv lig_charge    0
setenv withligand    1 
setenv myligand      'MyLigand.mol2' #this is the ligand that 

# choose from (0, 1]
# 0 = do not add water into binding site
# 1 = solvate occluded binding site (w/ 3D-RISM & GAsol)
setenv occluded_bs   1

setenv truncate      1
setenv truncate_dist 15
# truncate the system by 15 A around MyBindingSite.pdb

# length of MD simulation in nano second
setenv MD_ns         10


#=========================#
#    system preparation   #
#=========================#

if ( ! -f prot_amber.prmtop ) then

# step 1: protein residue protonation
if ( $protonation == 0 ) then
    python $WATSITEHOME/bin/run_reduce.py -p $myprotein
else if ( $protonation == 1 ) then
    /usr/bin/python $WATSITEHOME/bin/truncate_prot.py --pdb_propka $myprotein
else
    cp MyProtein.pdb prot.pdb
endif


# step 2: prepare ligand (if desired)
if ( $withligand == 1) then
    python $WATSITEHOME/bin/prepareLig_bcc.py -l $myligand -c $lig_charge
    # ligprep folder will be generated, lig_bcc_gaff.mol2 lig_AC.frcmod complex.pdb
endif


# step 3: solvate occluded binding site (w/ 3D-RISM & GAsol)
if ( $occluded_bs == 1) then

    if ( $withligand == 0) then
cat <<EOF > remove_wat.pml
load prot.pdb
remove solvent
remove (hydro)
save prot_dry.pdb, prot
quit
EOF
cat <<EOF > tleap_rism3d.in
source leaprc.protein.ff14SB
source leaprc.water.spce

prot = loadpdb prot_dry.pdb
saveamberparm prot prot_dry.prmtop prot_dry.inpcrd
savepdb prot prot_dry.pdb
quit
EOF

    else 
cat <<EOF > remove_wat.pml
load complex.pdb
remove solvent
remove (hydro)
save prot_dry.pdb, complex
quit
EOF
cat <<EOF > tleap_rism3d.in
source leaprc.protein.ff14SB
source leaprc.water.spce
source leaprc.gaff

SUB = loadmol2  ligprep/lig_bcc_gaff.mol2
loadamberparams ligprep/lig_AC.frcmod

prot = loadpdb prot_dry.pdb
saveamberparm prot prot_dry.prmtop prot_dry.inpcrd
savepdb prot prot_dry.pdb
quit
EOF

    endif
    
    pymol -c remove_wat.pml
    tleap -s -f tleap_rism3d.in
    bash  /programs/WATsite/bin/run_water_placement.sh

endif


# step 4: truncate system
if ( $truncate == 1) then

    if ( $occluded_bs == 1) then
        /usr/bin/python /programs/WATsite/bin/truncate_prot.py \
        --pdb_amber prot_dry.pdb --bs MyBindingSite.pdb -t $truncate_dist

    else
        if ( $withligand == 0) then
cat <<EOF > tleap_addH.in
source leaprc.protein.ff14SB
source leaprc.water.spce

prot = loadpdb prot.pdb
savepdb prot protH.pdb
quit
EOF
        else 
cat <<EOF > tleap_addH.in
source leaprc.protein.ff14SB
source leaprc.water.spce
source leaprc.gaff

SUB = loadmol2  ligprep/lig_bcc_gaff.mol2
loadamberparams ligprep/lig_AC.frcmod

prot = loadpdb complex.pdb
savepdb prot protH.pdb
quit
EOF
        endif
        tleap -s -f tleap_addH.in

        /usr/bin/python /programs/WATsite/bin/truncate_prot.py \
        --pdb_amber protH.pdb --bs MyBindingSite.pdb -t $truncate_dist

    endif
endif


# step 5: generate amber topology and coordinate endifles to run openmm simulation 

cat <<EOF > tleap_MD.in
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.spce
EOF

if ( $withligand == 1) then
cat <<EOF >> tleap_MD.in
SUB = loadmol2  ligprep/lig_bcc_gaff.mol2
loadamberparams ligprep/lig_AC.frcmod
EOF
endif        

if ( $truncate == 1 ) then
    if ( $occluded_bs == 1 ) then
cat <<EOF >> tleap_MD.in
dry = loadpdb prot_prep.pdb
wat  = loadpdb gasol_wat.pdb
prot = combine {dry wat}
EOF
    else
cat <<EOF >> tleap_MD.in
prot = loadpdb prot_prep.pdb
EOF
    endif
else
    if ( $occluded_bs == 1 ) then
cat <<EOF >> tleap_MD.in
dry = loadpdb prot_dry.pdb
wat  = loadpdb gasol_wat.pdb
prot = combine {dry wat}
EOF
    else
        if ( $withligand == 1) then
cat <<EOF >> tleap_MD.in
prot = loadpdb complex.pdb
EOF
        else
cat <<EOF >> tleap_MD.in
prot = loadpdb prot.pdb
EOF
        endif
    endif
endif

cat <<EOF >> tleap_MD.in
solvateBox prot SPCBOX 10.00
addIonsRand prot Cl- 0
addIonsRand prot Na+ 0
saveAmberParm prot prot_amber.prmtop prot_amber.inpcrd
quit

EOF
        
tleap -s -f tleap_MD.in 
ambpdb -p prot_amber.prmtop -c prot_amber.inpcrd > prot_amber.pdb

endif


#=========================#
#    run MD simulation    #
#=========================#

# prepare the system from amber endifles
if ( ! -f sys_min.xml ) then
/programs/miniconda/bin/python /programs/WATsite/bin/openmm_prep.py  -i ./prot_amber.inpcrd -t ./prot_amber.prmtop \
-l PME -c 10.00 --shake \
--restrain-mask '\!:WAT&\!@H=' -k 10.0 --reference prot_amber.pdb \
--cuda 0 
endif

# temperature coupling 298.15 K: 25,000 * 1.0 fs = 25 ps
if ( ! -f sys_NVT.xml ) then
/programs/miniconda/bin/python /programs/WATsite/bin/openmm_md.py  -i ./prot_amber.inpcrd -t ./prot_amber.prmtop \
--xml system.xml -s sys_min.xml --restart sys_NVT.xml \
-x sys_NVT.nc -r sys_NVT.info -o sys_NVT.out \
--temp 298.15 --gamma_ln 1.0 -n 25000 --interval 1000 --dt 1.0 \
--restrain-mask '\!:WAT&\!@H=' -k 2.5 --reference prot_amber.pdb \
--cuda 0 
endif

# pressure coupling 298.15 K and 1 bar: 1000000 (stpes) * 2.0 (timestep fs) = 2000.0 ps = 2.0 ns
if ( ! -f sys_NPT.xml ) then
/programs/miniconda/bin/python /programs/WATsite/bin/openmm_md.py  -i ./prot_amber.inpcrd -t ./prot_amber.prmtop \
--xml system.xml -s sys_NVT.xml --restart sys_NPT.xml \
-x sys_NPT.nc -r sys_NPT.info -o sys_NPT.out \
--temp 298.15 --gamma_ln 1.0 -n 1000000 --interval 10000 --dt 2.0 \
--restrain-mask '\!:WAT&\!@H=' -k 2.5 --reference prot_amber.pdb \
--cuda 0 --npt
endif

if ( ! -f sys_md_1.xml ) then
/programs/miniconda/bin/python /programs/WATsite/bin/openmm_md.py  -i ./prot_amber.inpcrd -t ./prot_amber.prmtop \
--xml system.xml -s sys_NPT.xml --restart sys_md_1.xml \
-x sys_md_1.nc  -r sys_md_1.info  -o sys_md_1.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '\!:WAT&\!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 0 -p prot_amber.pdb --cuda 0 --npt
endif

@ i = 2
while ($i <= $MD_ns)
    @ j = $i - 1
    @ k = ($i - 1) * 500

    if ( ! -f sys_md_${i}.xml ) then
    /programs/miniconda/bin/python /programs/WATsite/bin/openmm_md.py  -i ./prot_amber.inpcrd -t ./prot_amber.prmtop \
    --xml system.xml -s sys_md_${j}.xml --restart sys_md_${i}.xml \
    -x sys_md_${i}.nc  -r sys_md_${i}.info  -o sys_md_${i}.out \
    --temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
    --restrain-mask '\!:WAT&\!@H=' -k 2.5 --reference prot_amber.pdb \
    --water SPC/E --prev_frame $k -p prot_amber.pdb --cuda 0 --npt 
    endif

    @ i++ 
end



#================================#
#    hydration site prediction   #
#================================#

pymol -c /programs/WATsite/dat/super.pml
/programs/miniconda/bin/cpptraj < /programs/WATsite/dat/cpptraj.in

@ num_frame = $MD_ns * 500

cat <<EOF > input.bcf
\$Input_files:
   MyBindingSite.pdb          |input pdb file to define protein binding site    
      ../mdsnaps.pdb          |input trajectory in pdb                          
      ../waterEnergy          |folder contain water energy files                

\$Input_MDsnapshots:
               $num_frame          |number of frames in trajectory files             
               SPC/E          |water model (SPC/E; TIP3P; OPC; TIP4P; TIP4PEW; AMOEBA)

\$Grid_Parameters:
                 1.0          |size of water (SD)                               
                0.25          |griddelta: distance between two adjacent grid points
                 3.0          |Farest Distance to extend binding site box       
               0.045          |Water density for grid cutoff                    

\$Cluster_method:
                   2          |clustering method for HS identification (=1: DBSCAN; =2: QT clustering)

\$DBSCAN:
                 300          |min of points define neighbor in DBSCAN (start number)
                  80          |min of points define neighbor in DBSCAN (end number)

\$Hydro_Cluster:
                  10          |maximum clusters                                 
                2.50          |cluster_mean_dist                                
                2.75          |maxdist                                          
                10.0          |distance cutoff for QT clustering                

\$Entropy:
                   9          |covariance dimension (=3: two 3x3 matrix; =6: one 6x6 dimension; =9: both 3x3 and 6x6 method)
                  70          |number of bins used to construct Probability Distribution function

EOF

/programs/WATsite/bin/watsite  -c ./input.bcf -o WATsite_OUT

