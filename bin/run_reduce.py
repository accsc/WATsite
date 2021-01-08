import os, sys, re, getopt
WATsite_home = os.environ['WATSITEHOME']
reduce_exe_dir = WATsite_home + "/bin"
def run_reduce(argv):
	try:                                
		opts, args = getopt.getopt(argv, "-hp:r:") 
	except getopt.GetoptError:          
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-p':	
			inprot = arg
	#run reduce 
	comm1 = "%s/reduce -build %s > reduce_prep.pdb\n"%(reduce_exe_dir, inprot)
	print comm1
	os.system(comm1)
	
	# This is a check to see if the reduce_prep.pdb file is there
	try:
		Ri = open("reduce_prep.pdb", "r")
		try:
			pass
		finally:
			Ri.close()
	except IOError:
			print("\nError: can not find the file reduce_prep.pdb\n")
			sys.exit(0)
	
	# Open the reduce_prep.pdb file, obtain HIS protonation states, sort by res#, write out logfile
	Ri = open("reduce_prep.pdb", "r")
	
	ci = 0
	while Ri:
		line = Ri.readline()
		if line.find("USER  MOD") > -1:
			if line[25:28] == "HIS":
				ci += 1
		if line == '':
			break

	# Rewind, allocate the arrays to store res_num & protonation & order
	resnum = [0 for i in range(ci)]
	protna = [[] for i in range(ci)]

	
	Ri.seek(0)
	
	#store data to the arrays
	j = 0
	while Ri:
		line = Ri.readline()
		#for i in range(ci)
		if line.find("USER  MOD") > -1:
			if line[25:28] == "HIS":
				protonation = line[38:45].strip()
				resnum[j] = int(line[20:24].strip())
				if protonation == "no HD1":
					protna[j] = "HIE"
					j += 1
				elif protonation == "no HE2":
					protna[j] = "HID"
					j += 1
				elif protonation == "+bothHN":
					protna[j] = "HIP"
					j += 1
		if line == '':
			break
	
	# write prot.pdb for gromacs simulation
	Ri.seek(0)
	# Regular Expression (creates a pattern for the hydrogens)
	p = re.compile('H[A-Z]+')
	fpdb = open("prot.pdb", "w")
	
	m = -1
	cur_num = 0
	prev_num = -1
	for line in Ri:
		# if the line does not have a hydrogen
		if (line[0:4] == "ATOM" or line[0:6] == "HETATM") and p.match(line[12:16]) == None and p.match(line[13:16]) == None and line[13] != "H" and (line[13] != ' ' or (line[14:16] != 'H1' and line[14:16] != 'H2')):
			if line[17:20] == 'HIS':
				cur_num = int(line[22:26])
				if cur_num in resnum:
					if prev_num != cur_num:
						m += 1
						prev_num = cur_num
						line_1 = line[0:17] + protna[m] + " " + line[21:]
					else:
						line_1 = line[0:17] + protna[m] + " " + line[21:]
				else:
					line_1 = line
			else:
				line_1 = line
			fpdb.write(line_1)
		elif line[0:3] == "TER":
			fpdb.write(line)
		elif line[0:3] == "END":
			fpdb.write(line)
			break
	Ri.close()
	fpdb.close()
	os.remove("reduce_prep.pdb")

run_reduce(sys.argv[1:])
