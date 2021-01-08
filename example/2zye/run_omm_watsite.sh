# prepare the system from amber files
if [ ! -f sys_min.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_prep.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
-l PME -c 10.00 --shake \
--restrain-mask '!:WAT&!@H=' -k 10.0 --reference prot_amber.pdb \
--cuda 1 
fi

# temperature coupling 298.15 K: 25,000 * 1.0 fs = 25 ps
if [ ! -f sys_NVT.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_min.xml --restart sys_NVT.xml \
-x sys_NVT.nc -r sys_NVT.info -o sys_NVT.out \
--temp 298.15 --gamma_ln 1.0 -n 25000 --interval 1000 --dt 1.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--cuda 1 
fi

# pressure coupling 298.15 K and 1 bar: 1000000 (stpes) * 2.0 (timestep fs) = 2000.0 ps = 2.0 ns
if [ ! -f sys_NPT.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_NVT.xml --restart sys_NPT.xml \
-x sys_NPT.nc -r sys_NPT.info -o sys_NPT.out \
--temp 298.15 --gamma_ln 1.0 -n 1000000 --interval 10000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--cuda 1 --npt
fi

# production run 298.15 K and 1 bar: 10000000 (stpes) * 2.0 (timestep fs) = 20000.0 ps = 20.0 ns
if [ ! -f sys_md_1.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_NPT.xml --restart sys_md_1.xml \
-x sys_md_1.nc  -r sys_md_1.info  -o sys_md_1.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 0 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_2.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_1.xml --restart sys_md_2.xml \
-x sys_md_2.nc  -r sys_md_2.info  -o sys_md_2.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_3.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_2.xml --restart sys_md_3.xml \
-x sys_md_3.nc  -r sys_md_3.info  -o sys_md_3.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 1000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_4.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_3.xml --restart sys_md_4.xml \
-x sys_md_4.nc  -r sys_md_4.info  -o sys_md_4.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 1500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_5.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_4.xml --restart sys_md_5.xml \
-x sys_md_5.nc  -r sys_md_5.info  -o sys_md_5.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 2000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_6.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_5.xml --restart sys_md_6.xml \
-x sys_md_6.nc  -r sys_md_6.info  -o sys_md_6.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 2500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_7.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_6.xml --restart sys_md_7.xml \
-x sys_md_7.nc  -r sys_md_7.info  -o sys_md_7.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 3000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_8.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_7.xml --restart sys_md_8.xml \
-x sys_md_8.nc  -r sys_md_8.info  -o sys_md_8.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 3500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_9.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_8.xml --restart sys_md_9.xml \
-x sys_md_9.nc  -r sys_md_9.info  -o sys_md_9.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 4000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_10.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_9.xml --restart sys_md_10.xml \
-x sys_md_10.nc  -r sys_md_10.info  -o sys_md_10.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 4500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_11.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_10.xml --restart sys_md_11.xml \
-x sys_md_11.nc  -r sys_md_11.info  -o sys_md_11.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 5000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_12.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_11.xml --restart sys_md_12.xml \
-x sys_md_12.nc  -r sys_md_12.info  -o sys_md_12.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 5500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_13.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_12.xml --restart sys_md_13.xml \
-x sys_md_13.nc  -r sys_md_13.info  -o sys_md_13.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 6000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_14.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_13.xml --restart sys_md_14.xml \
-x sys_md_14.nc  -r sys_md_14.info  -o sys_md_14.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 6500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_15.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_14.xml --restart sys_md_15.xml \
-x sys_md_15.nc  -r sys_md_15.info  -o sys_md_15.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 7000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_16.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_15.xml --restart sys_md_16.xml \
-x sys_md_16.nc  -r sys_md_16.info  -o sys_md_16.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 7500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_17.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_16.xml --restart sys_md_17.xml \
-x sys_md_17.nc  -r sys_md_17.info  -o sys_md_17.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 8000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_18.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_17.xml --restart sys_md_18.xml \
-x sys_md_18.nc  -r sys_md_18.info  -o sys_md_18.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 8500 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_19.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_18.xml --restart sys_md_19.xml \
-x sys_md_19.nc  -r sys_md_19.info  -o sys_md_19.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 9000 -p prot_amber.pdb --cuda 1 --npt 
fi

if [ ! -f sys_md_20.xml ]; then
/programs/miniconda/bin//python /programs/WATsite/bin/openmm_md.py  -i /programs/WATsite/example/4nva_output/prot_amber.inpcrd -t /programs/WATsite/example/4nva_output/prot_amber.prmtop \
--xml system.xml -s sys_md_19.xml --restart sys_md_20.xml \
-x sys_md_20.nc  -r sys_md_20.info  -o sys_md_20.out \
--temp 298.15 --gamma_ln 1.0 -n 500000 --interval 1000 --dt 2.0 \
--restrain-mask '!:WAT&!@H=' -k 2.5 --reference prot_amber.pdb \
--water SPC/E --prev_frame 9500 -p prot_amber.pdb --cuda 1 --npt 
fi


