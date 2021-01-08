sed -i "s/HETATM/ATOM\ \ /g" MyBindingSite.pdb

center=`python /programs/WATsite/bin/get_centroid_pdb.py MyBindingSite.pdb`
echo $center

if [ ! -f prot_3drism.O.1.dx ]; then
/programs/rism3d/rism3d.snglpnt --pdb prot_dry.pdb --prmtop prot_dry.prmtop \
--xvv /programs/WATsite/dat/SPC_NaCl.xvv --guv prot_3drism \
--verbose 3 --volfmt dx --noasympcorr --maxstep 10000 --tolerance 1e-08 \
--ng 80,80,80 --solvbox 40,40,40 \
--centering 5 --centerM $center 
fi


IFS="," read -r -a xyz <<< "$center"
echo ${xyz[0]} 
echo ${xyz[1]} 
echo ${xyz[2]}

/programs/GAsol/gasol --dx prot_3drism.O.1.dx \
-r 10.0 -x ${xyz[0]} -y ${xyz[1]} -z ${xyz[2]} > gasol_wat.pdb

