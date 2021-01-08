#
# WATsite.py
#

import os,math,re
import string
import sys

if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import *
else:
    import tkinter as Tkinter
    from tkinter import *

import Pmw
#import distutils.spawn # used for find_executable
#import traceback

#from Tkinter import *
from tkFileDialog import *
from pymol import cmd, stored, util
from pymol.wizard import Wizard
from ScrolledText import *

import tkSimpleDialog
import tkMessageBox
import shutil
import subprocess
import colorsys
from collections import defaultdict
import numpy as np

monitor_list = []

global angstrom
global delta
angstrom = u"\u212B"
delta = u"\u0394" #unicode (somehow does not work on cheeta)
#delta = 'd'

global MAL_INFO 
MAL_INFO = "\nPyMOL WATsite Plugin: Hydration Site Prediction\nYing Yang, Amr Abdallah, Bingjie Hu, Markus A. Lill\nhttp://people.pharmacy.purdue.edu/~mlill/software/watsite/\n"



#######################################################################################
def __init__(self):
    cmd.set("retain_order", 1)
    cmd.set("pdb_use_ter_records", 1)

    # Simply add the menu entry and callback
    self.menuBar.addmenu('WATsite', 'Hydration Site Analysis',tearoff=TRUE)
    
    self.menuBar.addmenuitem('WATsite', 'command', 'Modify Paths',
        label = 'Modify Paths to Installed Programs',command = lambda s=self : modify_settings(s))
    self.menuBar.addmenuitem('WATsite', 'command', 'Prepare target',
        label = 'Prepare Protein/Ligand System',command = lambda s=self : prepare_protein(s))
    self.menuBar.addmenuitem('WATsite', 'command', 'Set Simulation Parameters',
        label = 'Set Parameters for Simulation',command = lambda s=self : setParam4MD(s))

    self.menuBar.addmenuitem('WATsite', 'separator')

    self.menuBar.addmenuitem('WATsite', 'command', 'Perform WATsite Analysis',
        label = 'Perform WATsite Analysis',command = lambda s=self : runWATsite(s))

    self.menuBar.addmenuitem('WATsite', 'separator')

    self.menuBar.addmenuitem('WATsite', 'command', 'Import/Monitor Results',
                label = 'Import WATsite Results',command = lambda s=self : import_results(s))
    self.menuBar.addmenuitem('WATsite', 'command', 'Estimate desolvation energy for ligand',
                label = 'Estimate Protein Desolvation Energy for Ligand',command = lambda s=self : estimate_desolvation(s))



#=================================================================================#
#=====================   Main: read from a previous setting   ====================#
#=================================================================================#
def read_settings():
    global python_dir
    global amber_home
    global pymol_dir
    global watsite_dir
    global reduce_dir

    """
    module load openmm-watsite
    prepend-path    PATH                   /usr/local/openmm-watsite/bin
    prepend-path    LD_LIBRARY_PATH        /usr/local/openmm-watsite/lib
    prepend-path    LD_LIBRARY_PATH        /usr/local/openmm-watsite/lib/plugins
    setenv          OPENMM_PLUGIN_DIR      /usr/local/openmm-watsite/lib/plugins
    setenv          OPENMM_CUDA_COMPILER   /usr/local/cuda-8.0/bin/nvcc

    module load anaconda
    prepend-path    PATH    /root/anaconda2/bin

    module load cuda8
    prepend-path    PATH                    /usr/local/cuda-8.0/bin
    prepend-path    LIBRARY_PATH            /usr/local/cuda-8.0/lib64
    prepend-path    LD_LIBRARY_PATH         /usr/local/cuda-8.0/lib64
    setenv          CUDA_HOME               /usr/local/cuda-8.0
    setenv          OPENMM_CUDA_COMPILER    /usr/local/cuda-8.0/bin/nvcc
    """
    
    path1 = os.environ.get('HOME')
    filename = "%s/WATsite_Settings.txt" % (path1)
    
    try:
        fi = open(filename, "r")
    except:
        tkMessageBox.showwarning("Settings file not found", "%s not found.\n You need to modify installed program paths first."% (filename) )
        #sys.exit()

    for i in fi:
        if i.find("python_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            python_dir = dat2.rstrip()
        if i.find("amber_home") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            amber_home = dat2.rstrip()
        if i.find("pymol_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            pymol_dir = dat2.rstrip()
        if i.find("watsite_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            watsite_dir = dat2.rstrip()
        if i.find("reduce_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            reduce_dir = dat2.rstrip()
    fi.close()
    os.environ["AMBERHOME"] = amber_home
    os.environ["PATH"] += os.pathsep + amber_home + '/bin'
    if not 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = ":"+amber_home + '/lib'
    elif not amber_home+'/lib' in os.environ.get("LD_LIBRARY_PATH"):
        os.environ['LD_LIBRARY_PATH'] += ":"+amber_home + '/lib'

    """
    print python_dir
    print amber_home
    print pymol_dir
    print watsite_dir
    print reduce_dir
    """


#===================== Menu 1: Modify paths to installed programs ====================#
class modifySettingsFile:
    def __init__(self, top):
        # Create the dialog.
        self.dialog = Pmw.Dialog(top,
                                 buttons = ('Exit','Save','Reset to defaults'),
                                 defaultbutton = 'Exit',
                                 title = 'PyMOL WATsite Plugin',
                                 command = self.apply)
        self.curdir = os.getcwd()
        
        w1 = Tkinter.Label(self.dialog.interior(),
                          text = MAL_INFO,
                          background = 'black',
                          foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='we')
        w1.grid_columnconfigure(0, weight=1)

        master = self.dialog.interior()

        Label(master, text='\nModify Paths to Installed Programs\n', foreground = 'blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=5)

        Tkinter.Label(master, text='python_dir\n(path to python associated with OpenMM installation):', anchor=W).grid(row=2, column=0)
        self.e_python_dir = Tkinter.Entry(master, width=70)
        self.e_python_dir.insert(0, python_dir)
        self.e_python_dir.grid(row=2, column=1)
        b_python_dir = Tkinter.Button(master, text='Browse', command = lambda: self.searchFOLDER(python_dir,self.e_python_dir)).grid(row=2, column=2)

        #Tkinter.Label(master, text='   ').grid(row=3, column=0)
        #Tkinter.Label(master, text='   ').grid(row=5, column=0)
        #Tkinter.Label(master, text='   ').grid(row=7, column=0)
        #Tkinter.Label(master, text='   ').grid(row=9, column=0)
        Tkinter.Label(master, text='   ').grid(row=20, column=0)

        Tkinter.Label(master, text='amber_home (main amber directory = $AMBERHOME):', anchor=W).grid(row=4, column=0)
        self.e_amber_home = Tkinter.Entry(master, width=70)
        self.e_amber_home.insert(0, amber_home)
        self.e_amber_home.grid(row=4, column=1)
        b_amber_home = Tkinter.Button(master, text='Browse', command = lambda: self.searchFOLDER(amber_home,self.e_amber_home)).grid(row=4, column=2)

        Tkinter.Label(master, text='pymol_dir (directory containing pymol executable):', anchor=W).grid(row=6, column=0)
        self.e_pymol_dir = Tkinter.Entry(master, width=70)
        self.e_pymol_dir.insert(0, pymol_dir)
        self.e_pymol_dir.grid(row=6, column=1)
        b_pymol_dir = Tkinter.Button(master, text='Browse', command = lambda: self.searchFOLDER(pymol_dir,self.e_pymol_dir)).grid(row=6, column=2)
        
        Tkinter.Label(master, text='watsite_dir (WATsite directory):', anchor=E).grid(row=8, column=0)
        self.e_watsite_dir = Entry(master, width=70)
        self.e_watsite_dir.insert(0, watsite_dir)
        self.e_watsite_dir.grid(row=8, column=1)
        self.b_watsite_dir = Button(master, text='Browse', command = lambda: self.searchFOLDER(watsite_dir,self.e_watsite_dir)).grid(row=8, column=2)

        #useReduce = IntVar()
        #Tkinter.Checkbutton(master, text="reduce_dir (directory containing reduce executable):", variable=useReduce, anchor=W).grid(row=6, column=0)
        Tkinter.Label(master, text='reduce_dir (directory containing reduce executable):', anchor=W).grid(row=10, column=0)
        self.e_reduce_dir = Tkinter.Entry(master, width=70)
        self.e_reduce_dir.insert(0, reduce_dir)
        self.e_reduce_dir.grid(row=10, column=1)
        b_reduce_dir = Tkinter.Button(master, text='Browse', command = lambda: self.searchFOLDER(reduce_dir,self.e_reduce_dir)).grid(row=10, column=2)

        self.dialog.activate(geometry = 'centerscreenalways')

    def reset(self, textcur, filecur):
        textcur.delete(0, END)
        textcur.insert(0, filecur)

    def writeNewFile(self):
        path_home = os.environ.get('HOME')
        filename = "%s/WATsite_Settings.txt" % path_home

        # write WATsite_Settings.txt file
        fi = open(filename, "w")
        fi.write("python_dir                     %s\n" % python_dir)
        fi.write("amber_home                     %s\n" % amber_home)
        fi.write("pymol_dir                      %s\n" % pymol_dir)
        fi.write("watsite_dir                    %s\n" % watsite_dir)
        fi.write("reduce_dir                     %s\n" % reduce_dir)
        fi.write("\n")
        fi.close()
        
    def searchFOLDER(self, filecur, textcur):
        self.folder_select = ""
        self.folder_select = askdirectory(title="Select directory", initialdir=filecur, mustexist=1)
        if self.folder_select:
            textcur.delete(0, END)
            textcur.insert(0, self.folder_select)
            filecur = self.folder_select

    def apply(self, result):
        global python_dir
        global amber_home
        global pymol_dir
        global watsite_dir
        global reduce_dir
        
        if result == 'Exit':
            self.dialog.deactivate()
            self.dialog.withdraw()
        elif result == 'Save':
            python_dir         = self.e_python_dir.get()
            amber_home         = self.e_amber_home.get()
            pymol_dir          = self.e_pymol_dir.get()
            watsite_dir        = self.e_watsite_dir.get()
            reduce_dir         = self.e_reduce_dir.get()
            self.writeNewFile()
            self.dialog.deactivate()
            self.dialog.withdraw()
        else:
            python_dir         = "/usr/local/anaconda2/bin/"
            amber_home         = "/usr/local/amber16"
            pymol_dir          = "/usr/local/pymol/"
            watsite_dir        = "/home/WATsite/"
            reduce_dir         = "/usr/local/reduce"
            self.reset(self.e_python_dir, python_dir)
            self.reset(self.e_amber_home, amber_home)
            self.reset(self.e_pymol_dir, pymol_dir)
            self.reset(self.e_watsite_dir, watsite_dir)
            self.reset(self.e_reduce_dir, reduce_dir)

def modify_settings(app):
    mfs = modifySettingsFile(app.root)


#=====================================================================================#
#=====================   Menu 2: Prepare Protien/Ligand System    ====================#
#=====================================================================================#
class prepareSystem4HydrationSite(tkSimpleDialog.Dialog):

    def __init__(self, top):
        # Create the dialog.
        self.dialog = Pmw.Dialog(top,
                                 buttons = ('Prepare System','Cancel'),
                                 defaultbutton = 'Cancel',
                                 title = 'PyMOL WATsite Plugin',
                                 command = self.apply)
        self.curdir = os.getcwd()
        
        w1 = Tkinter.Label(self.dialog.interior(),
                          text = MAL_INFO,
                          background = 'black', foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='wes')
        w1.grid_columnconfigure(0, weight=1)

        master = self.dialog.interior()
        curdir = os.getcwd()

        Label(master, text='\nSystem Preparation\n', foreground = 'blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=5)

        Tkinter.Label(master, text="Base directory").grid(row=2, column=0, sticky=W)
        self.eproj = Tkinter.Entry(master, width=45)
        self.eproj.insert(0, curdir)
        self.eproj.grid(row=2, column=1)
        self.bproj = Tkinter.Button(master, text='Browse', command=self.read_base_dir)
        self.bproj.grid(row=2, column=2)

        Tkinter.Label(master, text="Project subdirectory ").grid(row=3, column=0, sticky=W)
        self.eproj2 = Tkinter.Entry(master, width=45)
        self.eproj2.insert(0, "project_name")
        self.eproj2.grid(row=3, column=1)

        Tkinter.Label(master, text="  ").grid(row=4, columnspan=5, sticky=N)
        
        self.v0 = StringVar()
        self.rb_top1 = Tkinter.Radiobutton(master, text="Protein from file [PDB]", variable=self.v0, value="file_single")
        self.rb_top1.grid(row=5, column=0, sticky=W)
        self.e1_pdb = Tkinter.Entry(master, width=45)
        self.e1_pdb.insert(0, "Full path to protein file")
        self.e1_pdb.grid(row=5, column=1)
        self.e3_pdb = Tkinter.Button(master, text='Search and Import', command=self.searchPDB)
        self.e3_pdb.grid(row=5, column=2)
        
        self.rb_top3 = Tkinter.Radiobutton(master, text="Protein in current PyMol session", variable=self.v0, value="visual")
        self.rb_top3.grid(row=6, column=0, sticky=W)
        

        if len( cmd.get_names("objects") ) > 0:
            self.rb_top3.select()
        else:
            self.rb_top1.select()
        self.rb_sub0 = [[] for i in range(10000)]
        # top level radio box to select between ligand only, ligand + protein(fixed), ligand + protein(zone), all
        self.v_sub = StringVar()
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO')
            L = cmd.count_atoms('pro2')
            if L > 1:
                self.rb_sub0[ci] = Tkinter.Radiobutton(master, text=na, variable=self.v_sub, value=na)
                self.rb_sub0[ci].grid(row=ci+6, column=1, sticky=W)
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
        
        #self.rb_top3.select()
        if ci > 0:
            self.rb_sub0[0].select()

        #YY choice of reduce
        self.reduce = Tkinter.IntVar()
        Tkinter.Checkbutton(master, text="Use Reduce to optimize ASN,GLN,HIS conformations/protonation states", variable=self.reduce, onvalue=1, offvalue=0).grid(row=10, columnspan=5, sticky=W)
        Tkinter.Checkbutton(master, text="Rename residues to be AmberFF compatible based on hydrogens in PDB", variable=self.reduce, onvalue=2, offvalue=0).grid(row=11, columnspan=5, sticky=W)
        Tkinter.Checkbutton(master, text="Directly use AmberFF compatible residue names in PDB", variable=self.reduce, onvalue=3, offvalue=0).grid(row=12, columnspan=5, sticky=W)

        Tkinter.Label(master, text="  ").grid(row=20, columnspan=5, sticky=N)
        
        Tkinter.Label(master, text="Define binding site [PDB]").grid(row=21, column=0, sticky=W)
        self.bs_mol = Tkinter.Entry(master, width=45)
        self.bs_mol.insert(0, "Full path to binding site pdb file by known ligand or known binding site residues")
        self.bs_mol.grid(row=21, column=1)
        self.lig_mol = Button(master, text='Search and Import', command=self.searchBindingSiteFile)
        self.lig_mol.grid(row=21, column=2)
        Tkinter.Label(master, text="margin (%s)" % angstrom).grid(row=21, column=3)
        self.mb = Tkinter.Entry(master, width=5)
        self.mb.insert(0, "3.0")
        self.mb.grid(row=21, column=4)
                
        Tkinter.Label(master, text="  ").grid(row=22, columnspan=5, sticky=N)

        #YY choice of reduce
        self.occluded = Tkinter.IntVar()
        Tkinter.Checkbutton(master, text="Initial water placement with 3D-RISM and GAsol", variable=self.occluded, onvalue=1, offvalue=0).grid(row=23, columnspan=5, sticky=W)

        Tkinter.Label(master, text="  ").grid(row=24, columnspan=5, sticky=N)

        self.truncate = Tkinter.IntVar()
        Tkinter.Checkbutton(master, text="Truncate the protein from the center of the binding site", variable=self.truncate, onvalue=1, offvalue=0).grid(row=25, columnspan=2, sticky=W)
        Tkinter.Label(master, text="by cutoff (%s)" % angstrom).grid(row=25, column=3)
        self.truncate_cutoff = Tkinter.Entry(master, width=5)
        self.truncate_cutoff.insert(20, "20")
        self.truncate_cutoff.grid(row=25, column=4)

        Tkinter.Label(master, text="  ").grid(row=30, columnspan=5, sticky=N)

        self.wl = IntVar()
        self.wl_top3 = Tkinter.Radiobutton(master, text="Prediction without ligand", variable=self.wl, value=0)
        self.wl_top3.grid(row=31, column=0, sticky=W)               

        self.wl_top1 = Tkinter.Radiobutton(master, text="Prediction with ligand file [MOL2]", variable=self.wl, value=1)
        self.wl_top1.grid(row=32, column=0, sticky=W)
        self.wl_pdb = Tkinter.Entry(master, width=45)
        self.wl_pdb.insert(0, "Full path to ligand file")
        self.wl_pdb.grid(row=32, column=1)
        self.wl3_pdb = Tkinter.Button(master, text='Search and Import', command=self.searchLIGAND)
        self.wl3_pdb.grid(row=32, column=2)
        Tkinter.Label(master, text="ligand charge:").grid(row=32, column=3)
        self.lig_charge = Tkinter.Entry(master, width=5)
        self.lig_charge.insert(0, "0")
        self.lig_charge.grid(row=32, column=4)

        self.charge_method = Tkinter.IntVar()
        Tkinter.Checkbutton(master, text="Predict AM1-BCC charge charge using Antechamber", variable=self.charge_method, onvalue=1, offvalue=0).grid(row=33, column=1, sticky=W)     
        Tkinter.Checkbutton(master, text="Take user's charge in mol2 file", variable=self.charge_method, onvalue=2, offvalue=0).grid(row=34, column=1, sticky=W)

        Label(master, text="  ").grid(row=40, columnspan=5, sticky=N)
        self.ff_top1 = Tkinter.Label(master, text="Force Field")
        self.ff_top1.grid(row=41, column=0, sticky=W)
        self.ff = Tkinter.StringVar()
        self.ff.set('amber14SB')
        self.ff_top2 = Tkinter.OptionMenu(master, self.ff, 'amber14SB', 'amber99SBildn', 'amber99SB', 'AMOEBA-2013')
        self.ff_top2.grid(row=41, column=1, sticky=W)

        self.wat_top1 = Tkinter.Label(master, text="Water Model")
        self.wat_top1.grid(row=42, column=0, sticky=W)
        self.wat = Tkinter.StringVar()
        self.wat.set("SPC/E")
        #self.wat_top2 = Tkinter.OptionMenu(master, self.wat, 'SPC/E', 'TIP3P', 'OPC', 'TIP4P', 'TIP4PEW', 'AMOEBA03','AMOEBA14')
        self.wat_top2 = Tkinter.OptionMenu(master, self.wat, 'SPC/E', 'TIP3P', 'OPC', 'TIP4P', 'TIP4PEW', 'AMOEBA03')
        self.wat_top2.grid(row=42, column=1, sticky=W)
      
        self.box_top1 = Tkinter.Label(master, text="PBC Box Minimal Distance (%s)" % angstrom)
        self.box_top1.grid(row=43, column=0, sticky=W)
        self.box_top2 = Tkinter.Entry(master, width=9)
        self.box_top2.insert(0, "10.0")
        self.box_top2.grid(row=43, column=1, sticky=W)

        Tkinter.Label(master, text="  ").grid(row=50, columnspan=5, sticky=N)
        cmd.center()
        cmd.zoom()
        self.dialog.activate(geometry = 'centerscreenalways')

    def read_base_dir(self):
        base_dir = ""
        base_dir = askdirectory(title="Select base directory", mustexist=1)
        self.eproj.delete(0, END)
        self.eproj.insert(0, base_dir)   

    def searchPDB(self):
        #ftypes=(('pdb file', '*.pdb'), ('mol2 file', '*.mol2'), ('mol file', '*.mol'), ("All files", "*"))
        ftypes=(('pdb file', '*.pdb'), ("All files", "*"))
        indir = os.getcwd()
        self.g_pdb_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.g_pdb_name:
            self.e1_pdb.delete(0, END)
            self.e1_pdb.insert(0, self.g_pdb_name)
            cmd.load(self.g_pdb_name)

    def searchBindingSiteFile(self):
        ftypes=(('pdb file', '*.pdb'), ("All files", "*"))
        indir = os.getcwd()
        self.bs_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.bs_name:
            self.bs_mol.delete(0, END)
            self.bs_mol.insert(0, self.bs_name)
            cmd.load(self.bs_name)

    def searchLIGAND(self):
        ftypes=(('mol2 file', '*.mol2'), ("All files", "*"))
        indir = os.getcwd()
        self.ligand_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.ligand_name:
            self.wl_pdb.delete(0, END)
            self.wl_pdb.insert(0, self.ligand_name)
            cmd.load(self.ligand_name)
            #print self.ligand_name
        
    def searchFOLDER(self):
        self.folder_PDB = ""
        self.folder_PDB = askdirectory(title="Select directory", mustexist=1)
        if self.folder_PDB:
            self.e1_folder.delete(0, END)
            self.e1_folder.insert(0, self.folder_PDB)
    
    def getSystemFile(self):
        #if os.path.isdir(project_dir):
        #    print "Project directory exists. Delete and re-create one."
        #    os.system("rm -rf %s" % project_dir)
        if not os.path.isdir(project_dir):
            os.mkdir(project_dir)
        else:
            print "Project directory exists. Will overwrite previous files."

        os.chdir(project_dir)   

        if object_name == "visual":
            for na in cmd.get_names("objects"):
                if na == self.v_sub.get():
                    # remove hydrogens from protein
                    cmd.select('pro', na)
                    cmd.save("MyProtein.pdb", 'pro', 0, "pdb")
                    cmd.delete("pro")
        elif object_name == "file_single":
            os.system("cp %s ./MyProtein.pdb" % self.g_pdb_name)
        os.system("cp %s ./MyBindingSite.pdb" % self.bs_name)

    def runReduce(self):
        comm1 = "%s/reduce -build MyProtein.pdb > reduce_prep.pdb\n"%(reduce_dir)
        print "Running Commands: \n %s" % comm1
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
        Ri = open('reduce_prep.pdb', 'r')

        ci = 0
        cm = 0
        catoms = 0
        while Ri:
            line = Ri.readline()
            if line.find('USER  MOD S') > -1:
                #print line
                if line[25:28] == 'HIS':
                    ci += 1
                    #
                    flip_state = line[34:38].strip()
                    if flip_state == 'FLIP':
                        cm += 1
                elif line[25:28] == 'ASN':
                    flip_state = line[34:38].strip()
                    if flip_state == 'FLIP':
                        cm += 1
                elif line[25:28] == 'GLN':
                    flip_state = line[34:38].strip()
                    if flip_state == 'FLIP':
                        cm += 1
            catoms += 1
            if line == '':
                break
        # Rewind, allocate the arrays to store res_num & protonation & order
        order = [0 for i in range(ci)]
        resnum = [0 for i in range(ci)]
        protna = [[] for i in range(ci)]
        order2 = [0 for i in range(cm)]
        resnum2 = [0 for i in range(cm)]
        flipped = [[] for i in range(cm)]

        Ri.seek(0)
        # Store data to the arrays
        j = 0
        m = 0
        while Ri:
            line = Ri.readline()
            #for i in range (ci):
            if line.find('USER  MOD S') > -1:
                #print line
                if line[25:28] == 'HIS':
                    protonation =  line[38:45].strip()
                    resnum[j] = int(line[20:24].strip())
                    if protonation == 'no HD1':
                        protna[j]= 'HIE'
                        j += 1
                        flip_state = line[34:38].strip()
                        if flip_state == 'FLIP':
                            flipped[m] = 'HIE'
                            resnum2[m] = int(line[20:24].strip())
                            m += 1
                    elif protonation == 'no HE2':
                        protna[j]= 'HID'
                        j += 1 
                        flip_state = line[34:38].strip()
                        if flip_state == 'FLIP':
                            flipped[m] = 'HID'
                            resnum2[m] = int(line[20:24].strip())
                            m += 1
                    elif protonation == '+bothHN':
                        protna[j]= 'HIP' 
                        j += 1
                        flip_state = line[34:38].strip()
                        if flip_state == 'FLIP':
                            flipped[m] = 'HIP'
                            resnum2[m] = int(line[20:24].strip())
                            m += 1
                elif line[25:28] == 'ASN':
                    flip_state = line[34:38].strip()
                    if flip_state == 'FLIP':
                        flipped[m] = 'ASN'
                        resnum2[m] = int(line[20:24].strip())
                        m += 1
                elif line[25:28] == 'GLN':
                    flip_state = line[34:38].strip()
                    if flip_state == 'FLIP':
                        flipped[m] = 'GLN'
                        resnum2[m] = int(line[20:24].strip())
                        m += 1
            if line == '':
                break
        
        # Sort the array
        for i in range(ci):
            order[i] = i

        for i1 in range(ci):
            for i2 in range(i1+1, ci):
                if resnum[order[i2]] < resnum[order[i1]]:
                    tmpo = order[i1]
                    order[i1] = order[i2]
                    order[i2] = tmpo

        # Sort the array
        for i in range(cm):
            order2[i] = i

        for i1 in range(cm):
            for i2 in range(i1+1, cm):
                if resnum2[order2[i2]] < resnum2[order2[i1]]:
                    tmpo = order2[i1]
                    order2[i1] = order2[i2]
                    order2[i2] = tmpo

        # Open the log file and write 
        fo = open('Reduce.log', 'w')
        fo.write('Histidine states: \n')
        for i1 in range(ci):
            fo.write('%s %d\n' % (protna[order[i1]], resnum[order[i1]]))
        fo.write('\nHis, Asn or Gln sidechains flipped: \n')
        for i1 in range(cm):
            fo.write('%s %d\n' % (flipped[order2[i1]], resnum2[order2[i1]]))
        fo.close()
                       
        # Write out Latest .pdb 
        Ri.seek(0)
        # Regular Expression (creates a pattern for the hydrogens)
        p = re.compile('H[A-Z]+')
        fpdb = open('prot.pdb', 'w')

        cur_num = 0
        prev_num = -1
        m = -1
        for line in Ri:
            # If the line does not have a hydrogen
            if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') and p.match(line[12:16]) == None and p.match(line[13:16]) == None and line[13] != 'H' and (line[13] != ' ' or (line[14:16] != 'H1' and line[14:16] != 'H2')):
                if line[17:20] == 'HIS':
                    cur_num =  int(line[22:26])
                    if cur_num in resnum:
                        #print prev_num, cur_num
                        # Need a way to track residues.  Want to change the resdiue HIS to HIE not just a single line
                        if prev_num != cur_num:
                            m += 1
                            prev_num = cur_num
                            line_1 = line[0:17] + protna[order[m]] + ' ' + line[21:]
                        else:
                            line_1 = line[0:17] + protna[order[m]] + ' ' + line[21:]
                    else :
                            line_1 = line
                else :
                    line_1 = line
                fpdb.write(line_1)
        Ri.close()
        fpdb.close()
        os.remove("reduce_prep.pdb")

        Fi = open("prot.pdb", "r")
        Fo = open("tmp.pdb", "w")
        allLines = Fi.readlines()
        for i in range(len(allLines)):
                if allLines[i][0:6] == 'HETATM' and allLines[i-1][0:4] == 'ATOM':
                        Fo.write("TER\n")
                        Fo.write(allLines[i])
                else:
                        Fo.write(allLines[i])
        os.system("mv tmp.pdb prot.pdb")
 
    def rename_chargeRes(self):

        print "Rename residue based on hydrogens in input pdb"
        Fi = open("MyProtein.pdb", "r")
        chain_list = []
        for line in Fi:
            if line.find("ATOM") > -1 or line.find("HETATM") > -1: 
                if not line[21:22].upper() in chain_list:
                    chain_list.append(line[21:22].upper())
        print "Protein pdb (MyProtein.pdb) contains chain(s)", chain_list

        out = []
        for chain in chain_list:
            Fi.seek(0)
            print "Chain %s" % chain
            his_dict = defaultdict(list)
            asp_dict = defaultdict(list)
            glu_dict = defaultdict(list)
            lys_dict = defaultdict(list)

            for line in Fi:
                if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') and line[21:22].upper() == chain :
                    res_id = str( int(line[22:26]) )
                    atom_name = line[12:16].strip()
                    if line[17:20] == 'HIS':
                        his_dict[res_id].append(atom_name)
                    elif line[17:20] == 'ASP':
                        asp_dict[res_id].append(atom_name)
                    elif line[17:20] == 'GLU':
                        glu_dict[res_id].append(atom_name)
                    elif line[17:20] == 'LYS':
                        lys_dict[res_id].append(atom_name)
            #print his_dict, asp_dict, glu_dict, lys_dict
            protna = {}
            for key in his_dict:
                if "HD1" in his_dict[key] and "HE2" in his_dict[key]:
                    protna[key] = "HIP"
                elif "HD1" in his_dict[key] and not "HE2" in his_dict[key]:
                    protna[key] = "HID"
                elif "HE2" in his_dict[key] and not "HD1" in his_dict[key]:
                    protna[key] = "HIE"
                else:
                    protna[key] = "HIM"
            for key in asp_dict:
                if "HD2" in asp_dict[key]:
                    protna[key] = "ASH"
                else:
                    protna[key] = "ASP"
            for key in glu_dict:
                if "HE2" in glu_dict[key]:
                    protna[key] = "GLH"
                else:
                    protna[key] = "GLU"
            for key in lys_dict:
                if "HZ3" in lys_dict[key]:
                    protna[key] = "LYS"
                else:
                    protna[key] = "LYN"
            print protna

            Fi.seek(0)
            for line in Fi:
                if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') and line[21:22].upper() == chain :
                    res_id = str(int(line[22:26]))
                    atom_name = line[12:16].strip()
                    alternative_loc = str(line[16])
                    if atom_name[0] != "H" and alternative_loc != "B":
                        if res_id in protna:
                            out.append( line[0:17] + protna[res_id] + " " + line[21:] )
                        else:
                            out.append( line )
            out.append("TER\n")

        print 'HIS/ ASP/ GLU/ LYS residues are renamed'
        Fi.close()

        Fo = open("prot.pdb", "w")
        for i in range(len(out)):
            if out[i][0:6] == 'HETATM' and out[i-1][0:4] == 'ATOM':
                Fo.write("TER\n")
                Fo.write(out[i])
            elif out[i].find("HOH") > -1 and not out[i-1].find("HOH")>-1:
                Fo.write("TER\n")
                Fo.write(out[i])
            else:
                Fo.write(out[i])
        Fo.close()

    def runAcpype(self):
        curdir = os.getcwd()
        # first rename ligand res name to SUB
        if not os.path.isdir("ligPrep"):
            os.mkdir("ligPrep")
        os.chdir("ligPrep")

        fi = open(self.ligand_name, "r")
        fo = open("lig.mol2", "w")

        line = fi.readline()
        while line != '':
            if line.find("@<TRIPOS>MOLECULE") > -1:
                fo.write(line)
                line = fi.readline()
                fo.write(line)
                line = fi.readline()
                fo.write(line)
                numatoms = int(line.split()[0])
            elif line.find("@<TRIPOS>ATOM") > -1:
                fo.write(line)
                hcount = 1
                for i in range(numatoms):
                    line = fi.readline()
                    ss = line.split()
                    #if ss[1][0] == 'H':
                    #    atom_name = 'H%d'%hcount
                    #    hcount += 1
                    #else:
                    #    atom_name = ss[1]
                    atom_name = ss[1]
                    newline = "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n"%(int(ss[0]), atom_name, float(ss[2]), float(ss[3]), float(ss[4]), ss[5], float(ss[-1]))
                    fo.write(newline)
            else:
                fo.write(line)
            line = fi.readline()
        fi.close()
        fo.close()      

        #AmaberTools 
        """To call Antechamber and execute it

        Usage: antechamber -i  input file name
                           -fi input file format
                           -o  output file name
                           -fo output file format
                           -c  charge method
                           -cf charge file name
                           -nc net molecular charge (int)
                           -a  additional file name
                           -fa additional file format
                           -ao additional file operation
                                crd  only read in coordinate
                                crg only read in charge
                                name   only read in atom name
                                type   only read in atom type
                                bond   only read in bond type
                           -m  multiplicity (2S+1), default is 1
                           -rn residue name, overrides input file, default is MOL
                           -rf residue toplogy file name in prep input file,
                                      default is molecule.res
                           -ch check file name for gaussian, default is molecule
                           -ek mopac or sqm keyword, inside quotes
                           -gk gaussian keyword, inside quotes
                           -df am1-bcc flag, 2 - use sqm(default); 0 - use mopac
                           -at atom type, can be gaff (default), amber, bcc and sybyl
                           -du fix duplicate atom names: yes(y)[default] or no(n)
                           -j  atom type and bond type prediction index, default is 4
                                0     no assignment
                                1     atom type
                                2     full  bond types
                                3     part  bond types
                                4     atom and full bond type
                                5     atom and part bond type
                           -s  status information: 0(brief), 1(default) or 2(verbose)
                           -pf remove intermediate files: yes(y) or no(n)[default]
                           -i -o -fi and -fo must appear; others are optional

         file format type  abbre. index | file format type abbre. index
        ---------------------------------------------------------------
        Antechamber        ac       1  | Sybyl Mol2         mol2    2
        PDB                pdb      3  | Modified PDB       mpdb    4
        AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6
        Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8
        Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10
        Gaussian Output    gout    11  | Mopac Output       mopout 12
        Alchemy            alc     13  | CSD                csd    14
        MDL                mdl     15  | Hyper              hin    16
        AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18
        Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20
        Divcon Input       divcrt  21  | Divcon Output      divout 22
        Charmm             charmm  23  | SQM Output         sqmout 24
        Charmm             charmm  25  | Gaussian ESP       gesp   26
        Component cif      ccif    27  |                             
        --------------------------------------------------------------

                      List of the Charge Methods

        charge method     abbre.  index | charge method    abbre. index
        ----------------------------------------------------------------
        RESP               resp     1  |  AM1-BCC            bcc     2
        CM1                cm1      3  |  CM2                cm2     4
        ESP (Kollman)      esp      5  |  Mulliken           mul     6
        Gasteiger          gas      7  |  Read in charge     rc      8
        Write out charge   wc       9  |  Delete Charge      dc     10
        ----------------------------------------------------------------
        """
        #
        if os.path.isfile("lig_bcc_gaff.mol2"):
            print "\n\n"+"#"+"="*25+"#"
            print "Found output from previous antechamber run.\nDelete the previous project or create a new project if you want to perform antechamber again.\n"
        else:
            if charge_method == 2:
                cm = ''
            elif charge_method == 1:
                cm = '-c bcc'
            command = '%s/bin/antechamber -i lig.mol2 -fi mol2 -o lig_bcc_gaff.mol2 -fo mol2 %s -nc %i -m 1 -s 2 -df 2 -at gaff -pf y -j 1' % (amber_home, cm, lig_charge)
            print "Executing Antechamber..." 
            print command
            os.system(command)

        if not os.path.isfile("lig_bcc_gaff.mol2"):
            print "\nError: antechamber aborting ! "
            hint1 = "HINT1: is path of amber_home correctly set?\n"
            hint2 = "HINT2: is ligand charge correctly set?\n"
            print hint1,hint2
            sys.exit(1)

        if not os.path.isfile("lig_AC.frcmod"):
            command = '%s/bin/parmchk -i lig_bcc_gaff.mol2 -f mol2 -o lig_AC.frcmod' % (amber_home)
            print "Executing parmchk...\n %s" % command
            os.system(command)

        os.chdir(curdir)

        cmd.delete("all")
        cmd.load("prot.pdb")
        cmd.load("ligPrep/lig.mol2")
        cmd.save("complex.pdb", "prot or lig", 0, "pdb")
        cmd.delete("prot")
        cmd.delete("lig")
        cmd.load("complex.pdb")

        fo = open('tmp.pdb', 'w')
        fi = open('complex.pdb', 'r')
        lines = fi.readlines()
        fi.close()
        for i in range( len(lines) ):
            if lines[i].find("OXT") > -1:
                fo.write(lines[i])
                fo.write("TER\n")
            elif lines[i].find("HOH") > -1 and not lines[i-1].find("HOH") > -1:
                fo.write("TER\n")
                fo.write(lines[i])
            elif lines[i].find("SUB") > -1 and not lines[i-1].find("SUB") > -1:
                fo.write("TER\n")
                fo.write(lines[i])
            else:
                fo.write(lines[i])
        fo.close()
        os.system("mv tmp.pdb complex.pdb")

    def runTleapForAmoeba(self):
        print "\n"+"#"+"="*25+"#"
        print "Now use tleap to add hydrogens and solvate the system..."
        print "Use Force Field: %s" % forcefield
        print "Use Water Model: %s" % watermodel
        print "We will prepare the system with amber14SB force field and SPC/E water model, and convert to AMOEBA-2013 atom types with OpenMM."
        print "Create PBC Box by Extending %.2f angstrom" % (min_dis)
        print "#"+"="*25+"#"

        fo = open('tmp.pdb', 'w')
        fi = open('prot.pdb', 'r')
        for line in fi:
            if line.find("OXT") > -1:
                fo.write(line)
                fo.write("TER\n")
            else:
                fo.write(line)
        fo.close()
        fi.close()
        os.system("mv tmp.pdb prot.pdb")

        fp = open("tleap.in", "w")
        fp.write("source leaprc.protein.ff14SB\n")
        fp.write("source leaprc.gaff\n" )
        fp.write("source leaprc.water.spce\n")
        fp.write("prot = loadpdb prot.pdb\n")
        fp.write("solvateBox prot SPCBOX %.2f \n" % (min_dis) )
        fp.write("charge prot\n")
        fp.write("addIons2 prot Cl- 0\n")
        fp.write("addIons2 prot Na+ 0\n")
        fp.write("saveAmberParm prot prot_amber.prmtop prot_amber.inpcrd\n")
        fp.write("quit\n")
        fp.write("\n")
        fp.close()

        command = '%s/bin/tleap -f tleap.in' % (amber_home)
        print "Executing tleap...\n %s" % command
        os.system(command)

        command = '%s/bin/ambpdb -p prot_amber.prmtop -c prot_amber.inpcrd -ep -conect > prot_amber.pdb' % (amber_home)
        print "Executing ambpdb...\n %s" % command
        os.system(command)

        cmd.delete("all")
        cmd.load("prot_amber.pdb")

    def run_water_placement(self):
        fp = open("tleap_rism3d.in", "w")
        fp.write("source leaprc.protein.ff14SB\n")
        fp.write("prot = loadpdb prot_dry.pdb\n")
        fp.write("saveamberparm prot prot_dry.prmtop prot_dry.inpcrd\n")
        fp.write("savepdb prot prot_dry.pdb\n")
        fp.write("quit\n")
        fp.close()
        command = '%s/bin/tleap -f tleap_rism3d.in' % (amber_home)
        print "Executing...\n %s" % command
        os.system(command)

        command = 'bash %s/bin/run_water_placement.sh' % (watsite_dir)
        print "Executing...\n %s" % command
        subprocess.call(command, shell=True)

    def truncate_and_cap(self):
        with open("tleap_addH.in", "w") as fo:
            fo.write("source leaprc.protein.ff14SB\n")
            fo.write("source leaprc.water.spce\n")
            fo.write("prot = loadpdb prot.pdb\n")
            fo.write("savePdb prot protH.pdb\nquit\n")
        command = 'tleap -f tleap_addH.in'
        print "Executing...\n %s" % command
        os.system(command)

        lines = []
        chain = 0
        chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        with open("protH.pdb") as protH:
            for line in protH:
                if line.split()[0] == "TER":
                    chain += 1
                elif line.split()[0] == "ATOM" and line.split()[3] != "WAT":
                    line = line[:21] + chains[chain] + line[22:]
                lines.append(line)

        with open("protH.pdb", "w") as out:
            for line in lines:
                out.write(line)

        cmd.delete("all")
        cmd.load("protH.pdb")
        cmd.load("MyBindingSite.pdb")

        cmd.remove('(hydro)')

        cmd.select('br. MyBindingSite around %f and protH' % truncate_cutoff)
        cmd.save('prot_truncated.pdb', 'sele')
        cmd.load('prot_truncated.pdb')

        inpdb = 'prot_truncated.pdb'
        fullpdb = 'protH.pdb'
        outpdb = "prot_prep.pdb"

        info  = {}
        lines = defaultdict(list)
        with open(fullpdb, "r") as f:
            for line in f:
                if line.find("ATOM") > -1 or line.find("HETATM") > -1:
                    resNum, ele, x, y, z = int(line[22:26]), line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54])
                    if ele in ['C','CA','O','N']:
                        info['%d-%s' %(resNum, ele)] = np.array([x, y, z])
                    lines[resNum].append(line)

        fi = open(inpdb, "r")
        fo = open(outpdb, "w")
        resNumList = []
        for line in fi:
            if (line.find("ATOM") > -1 or line.find("HETATM") > -1) and not (line.find("WAT") > -1):
                resNum, ele = int(line[22:26]), line[12:16].strip()
                if resNum not in resNumList: # N atom of a new residue
                    resNumList.append(resNum)
                    if len(resNumList) == 1: # first residue
                        find_CA = "%d-CA" % (resNum-1)
                        find_C  = "%d-C"  % (resNum-1)
                        find_O  = "%d-O"  % (resNum-1)
                        if find_CA in info:
                            chain_id = lines[ resNum - 1 ][0][21:22]

                            fo.write("ATOM  %5d  CH3 ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) )
                            fo.write("ATOM  %5d  C   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_C][0], info[find_C][1], info[find_C][2]) )
                            fo.write("ATOM  %5d  O   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_O][0], info[find_O][1], info[find_O][2]) )
                        else:
                            ### MAY CAUSE PROBLEM
                            pass
                            # don't write at the beginning of original terminal # YY 3-27
                            #fo.write("ATOM  %5d  C   ACE  %4d    %8.3f%8.3f%8.3f\n" % (0, resNumList[-1], 0, 0, 0) )
                        fo.write(line)
                    else:
                        if resNumList[-1] == resNumList[-2] + 2: # only one residue deleted, add it back --> we don't break 
                            res_add = resNumList[-2] + 1
                            for line in lines[res_add]:
                                fo.write(line)
                            fo.write(line)
                        elif resNumList[-1] == resNumList[-2] + 1: # no break
                            fo.write(line)
                        else: # this resnum is larger than prev +1 
                            find_N  = "%d-N"  % (resNumList[-2] + 1)
                            find_CA = "%d-CA" % (resNumList[-2] + 1)
                            chain_id = lines[ resNumList[-2] + 1 ][0][21:22]
                            fo.write("ATOM  %5d    N NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-2], info[find_N ][0], info[find_N ][1], info[find_N ][2]) )
                            fo.write("ATOM  %5d  CH3 NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-2], info[find_CA][0], info[find_CA][1], info[find_CA][2]) ) 

                            fo.write("TER\n")

                            find_CA = "%d-CA" % (resNumList[-1]-1)
                            find_C  = "%d-C"  % (resNumList[-1]-1)
                            find_O  = "%d-O"  % (resNumList[-1]-1)
                            chain_id = lines[ resNumList[-1]-1 ][0][21:22]
                            fo.write("ATOM  %5d  CH3 ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) )
                            fo.write("ATOM  %5d  C   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_C ][0], info[find_C ][1], info[find_C ][2]) )
                            fo.write("ATOM  %5d  O   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_O ][0], info[find_O ][1], info[find_O ][2]) )

                            fo.write(line)
                else: # for other atoms in the same residue
                    fo.write(line)
            elif line.find("END") > -1:
                resNum = resNumList[-1] # prev res
                find_N  = "%d-N"  % (resNum+1)
                find_CA = "%d-CA" % (resNum+1)
                if find_CA in info:
                    chain_id = lines[ resNum+1 ][0][21:22]
                    fo.write("ATOM  %5d    N NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_N ][0], info[find_N ][1], info[find_N ][2]) )
                    fo.write("ATOM  %5d  CH3 NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) ) 
                else:
                    pass
                fo.write("END\n")
            elif line.find("TER") > -1:
                pass
            else: # not ATOM/HETATM
                fo.write(line)
        fi.close()
        fo.close()

        fo = open('tmp', "w")
        ter = 0
        with open(outpdb, "r") as fi:
            for line in fi:
                #if (line.find("OXT") > -1 or line.find(" H2 ") > -1 or line.find(" H3 ") > -1 ) and not line.find("WAT") > -1:
                # OXT or H1-3 should not be added if tleap runs correctly. // YY 3-28
                # It was added because not chain was specified in the cap.pdb file
                """
                if not line.find("WAT") > -1:
                    pass
                else:
                    if line.find("WAT") > -1 and ter==0:
                        fo.write("TER\n")
                        fo.write(line)
                        ter = 1
                    else:
                        fo.write(line)
                """
                if line.find("WAT") > -1 and ter==0:
                    fo.write("TER\n")
                    fo.write(line)
                    ter = 1
                else:
                    fo.write(line)

        fo.close()
        os.system("mv tmp %s" % outpdb)

    def runTleap(self):
        print "\n\n"+"#"+"="*25+"#"
        if occluded:
            print "3D-RISM and GAsol are used to place water molecules"
        if truncate:
            print "Truncate the protein from the center of the binding site by %.2f angstrom" % (truncate_cutoff)
        print "\n\n"+"#"+"="*25+"#"
        print "Now use tleap to add hydrogens and solvate the system..."
        print "Use Force Field: %s" % forcefield
        print "Use Water Model: %s" % watermodel
        print "Create PBC Box by Extending %.2f angstrom" % (min_dis)
        print "#"+"="*25+"#"

        fp = open("tleap_MD.in", "w")
        if forcefield == "amber14SB":
            fp.write("source leaprc.protein.ff14SB\n")
        else:
            fp.write("source oldff/leaprc.ff%s\n" % forcefield.split('amber')[-1] )
        fp.write("source leaprc.gaff\n" )
        if watermodel == "TIP4P":
            fp.write("loadOff atomic_ions.lib\n")
            fp.write("loadOff solvents.lib\n")
            fp.write("HOH = TP4\n")
            fp.write("WAT = TP4\n")
            fp.write("loadAmberParams frcmod.tip4p\n")
            fp.write("loadAmberParams frcmod.ionsjc_tip4pew\n")
            fp.write("loadAmberParams frcmod.ions234lm_126_tip4pew\n")
        elif watermodel == "SPC/E":
            fp.write("source leaprc.water.spce\n")
        else:
            fp.write("source leaprc.water.%s\n" % str.lower(watermodel) )
        if withligand:
            fp.write("loadamberparams ligPrep/lig_AC.frcmod\n")
            fp.write("SUB = loadmol2  ligPrep/lig_bcc_gaff.mol2\n")
            fp.write("prot = loadpdb  complex.pdb\n")
        else:
            if occluded and not truncate:
                fp.write("dry = loadpdb prot_dry.pdb\n")
                fp.write("wat  = loadpdb gasol_wat.pdb\n")
                fp.write("prot = combine {dry wat}\n")
            elif occluded and truncate:
                fp.write("dry = loadpdb prot_prep.pdb\n")
                fp.write("wat  = loadpdb gasol_wat.pdb\n")
                fp.write("prot = combine {dry wat}\n")
            elif not occluded and truncate:
                fp.write("prot = loadpdb prot_prep.pdb\n")
            else:
                fp.write("prot = loadpdb prot.pdb\n")

        if watermodel == "SPC/E":
            fp.write("solvateBox prot SPCBOX %.2f\n" % (min_dis) )
        else:
            fp.write("solvateBox prot %sBOX %.2f\n" % (watermodel, min_dis) )
        fp.write("addIons2 prot Cl- 0\n")
        fp.write("addIons2 prot Na+ 0\n")
        fp.write("saveAmberParm prot prot_amber.prmtop prot_amber.inpcrd\n")
        fp.write("quit\n")
        fp.write("\n")
        fp.close()

        if os.path.isfile("prot_amber.prmtop"):
            print "\n\n"+"#"+"="*25+"#"
            print "Found output from previous tleap run.\nDelete the previous project or create a new project if you want to perform tleap again.\n"
        else:
            command = '%s/bin/tleap -f tleap_MD.in' % (amber_home)

            print "Executing...\n %s" % command
            os.system(command)

            command = '%s/bin/ambpdb -p prot_amber.prmtop -c prot_amber.inpcrd -ep > prot_amber.pdb' % (amber_home)
            print "Executing...\n %s" % command
            os.system(command)

        cmd.load("prot_amber.pdb")

    def apply(self, result):
        global project_dir
        global project_dir_short
        global object_name

        global useReduce 
        global occluded 
        global truncate 
        global truncate_cutoff

        global useProtna
        global withligand
        global lig_charge
        global watermodel
        global forcefield
        global min_dis
        global margin
        global charge_method

        if result == 'Cancel':
            #self.parent.focus_set()
            self.dialog.deactivate()
            self.dialog.withdraw()
        elif result == 'Prepare System':
            #print "Now preparing system..."
            project_dir = ""
            project_dir_short = ""
            project_dir = "%s/%s" % (self.eproj.get(), self.eproj2.get())
            project_dir_short = "%s" % (self.eproj2.get())
            object_name = self.v0.get()


            useReduce = int(self.reduce.get())  # 1 (reduce) or 2 (propka hydrogens) or 3 (use residue name in PDB)
            occluded = int(self.occluded.get())  # 1 (truncate) or 0 (no truncation) 
            truncate = int(self.truncate.get())  # 1 (truncate) or 0 (no truncation) 
            truncate_cutoff = float(self.truncate_cutoff.get())

            #useProtna = self.protna.get()  # 0 or 1
            withligand = self.wl.get()     # 0 or 1
            lig_charge = int(self.lig_charge.get())
            charge_method = self.charge_method.get() 
            watermodel = self.wat.get()
            forcefield = self.ff.get()
            min_dis    = float(self.box_top2.get())
            margin     = float(self.mb.get())
            self.margin = float(self.mb.get())

            print "Creating the project directory: %s ..." % project_dir
            self.getSystemFile()

            fo = open("sys_prep", "w")
            fo.write("%f\n"%margin)
            fo.write("%s\n"%forcefield)
            fo.write("%s\n"%watermodel)          
            fo.close()

            if useReduce == 1:
                self.runReduce()
            elif useReduce == 2:  #if useProtna:
                self.rename_chargeRes()
            elif useReduce == 3: # directly use residue names in PDB
                print "Directly use pdb"
                os.system("cp MyProtein.pdb prot.pdb")
            else:
                print "Error: Please select how to treat protonation states..."
                sys.exit(-1)

            if withligand and forcefield != "AMOEBA-2013":
                self.runAcpype()
                os.system("cp %s ./MyLigand.mol2" % self.wl_pdb.get() )
                

            if occluded:
                cmd.delete("all")
                cmd.load("prot.pdb")
                cmd.remove("solvent")
                cmd.save("prot_dry.pdb", "prot", 0, "pdb")
                self.run_water_placement() # --> generated gasol_wat.pdb
                cmd.load("prot_dry.pdb")
                cmd.load("gasol_wat.pdb")

            if truncate:
                self.truncate_and_cap() # --> prot_prep.pdb


            self.runTleap()

            #if withligand and forcefield == "AMOEBA-2013": 
            #    tkMessageBox.showwarning("Warning","Prediction with ligand in AMOEBA force field is currently not supported...\n")
            #    #self.dialog.deactivate()
            #    self.dialog.withdraw()
            #    return -1
            #elif withligand and forcefield != "AMOEBA-2013":
            #    self.runAcpype()
            #    os.system("cp %s ./MyLigand.mol2" % self.wl_pdb.get() )

            #if forcefield == "AMOEBA-2013": 
            #    self.runTleapForAmoeba()
            #else:
            #    self.runTleap()

            self.dialog.deactivate()
            self.dialog.withdraw()

def prepare_protein(app):       
    prep = prepareSystem4HydrationSite(app.root)


#====================================================================================#
#=====================   Menu 3: Set Parameters for Simulation   ====================#
#====================================================================================#
class setSimulationParameter:
    def __init__(self, top):
        # Create the dialog.
        self.dialog = Pmw.Dialog(top,
                                 buttons = ('Save Parameters','Cancel'),
                                 defaultbutton = 'Cancel',
                                 title = 'PyMOL WATsite Plugin',
                                 command = self.apply)
        self.curdir = os.getcwd()
        
        w1 = Tkinter.Label(self.dialog.interior(),
                          text = MAL_INFO, 
                          background = 'black',
                          foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='we')
        w1.grid_columnconfigure(0, weight=1)

        master = self.dialog.interior()

        Label(master, text='\nSet Parameters for Simulation\n', foreground = 'blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=5)

        curdir = os.getcwd()
        Tkinter.Label(master, text="Project Directory").grid(row=2, column=0, sticky=W)
        self.proj_dir = Tkinter.Entry(master, width=50)
        try:
            project_dir
        except:
            self.proj_dir.insert(0, curdir)
        else:
            self.proj_dir.insert(0, project_dir)
        self.proj_dir.grid(row=2, column=1)
        self.bproj = Tkinter.Button(master, text='Browse', command=self.read_base_dir)
        self.bproj.grid(row=2, column=2)


        Tkinter.Label(master, text="Prepared Amber Topology File").grid(row=3, column=0, sticky=W)
        self.ambTop = Tkinter.Entry(master, width=50)
        try:
            project_dir
        except:
            self.ambTop.insert(0, "%s/prot_amber.prmtop" % curdir )
        else:
            self.ambTop.insert(0, "%s/prot_amber.prmtop" % project_dir )
        self.ambTop.grid(row=3, column=1)
        self.butTop = Tkinter.Button(master, text='Browse', command=self.searchTopFile )
        self.butTop.grid(row=3, column=2)

        Tkinter.Label(master, text="Prepared Amber Coordinate File").grid(row=4, column=0, sticky=W)
        self.ambCrd = Tkinter.Entry(master, width=50)
        try:
            project_dir
        except:
            self.ambCrd.insert(0, "%s/prot_amber.inpcrd" % curdir )
        else:
            self.ambCrd.insert(0, "%s/prot_amber.inpcrd" % project_dir )
        self.ambCrd.grid(row=4, column=1)
        self.butCrd = Tkinter.Button(master, text='Browse', command=self.searchCrdFile )
        self.butCrd.grid(row=4, column=2)

        Tkinter.Label(master, text="  ").grid(row=5, columnspan=5, sticky=N)

        self.cuda_top = Tkinter.Label(master, text="CUDA Platform DeviceIndex")
        self.cuda_top.grid(row=6, column=0, sticky=W)
        self.cuda_index = Tkinter.Entry(master, width=9)
        self.cuda_index.insert(0, "0")
        self.cuda_index.grid(row=6, column=1, sticky=W)

        Label(master, text="  ").grid(row=7, columnspan=5, sticky=N)

        self.lr_top1 = Tkinter.Label(master, text="Long-range Electrostatic Treatment")
        self.lr_top1.grid(row=8, column=0, sticky=W)
        self.lr = Tkinter.StringVar()
        self.lr.set('PME')
        self.lr_top2 = Tkinter.OptionMenu(master, self.lr, 'PME', 'Ewald', 'NoCutoff')
        self.lr_top2.grid(row=8, column=1, sticky=W)

        self.lr_top3 = Tkinter.Label(master, text="Non-bonded Interactions Cutoff (%s)" % angstrom)
        self.lr_top3.grid(row=9, column=0, sticky=W)
        self.lr_cut = Tkinter.Entry(master, width=9)
        self.lr_cut.insert(0, "10.0")
        self.lr_cut.grid(row=9, column=1, sticky=W)

        Label(master, text="  ").grid(row=10, columnspan=5, sticky=N)

        self.shake = Tkinter.IntVar()
        self.shake_top1 = Tkinter.Checkbutton(master, text="Apply constrain on bonds containing hydrogen (using SETTLE on waters)", variable=self.shake, onvalue=1, offvalue=0)
        self.shake.set(1)
        self.shake_top1.grid(row=11, columnspan=5, sticky=W)

        Label(master, text="  ").grid(row=12, columnspan=5, sticky=N)


        self.restrain = Tkinter.IntVar()
        self.restrain_top1 = Tkinter.Checkbutton(master, text="Apply positional restrain on non-hydrogen atoms", variable=self.restrain, onvalue=1, offvalue=0)
        self.restrain.set(1)
        self.restrain_top1.grid(row=13, columnspan=5, sticky=W)

        self.restrain_top2 = Tkinter.Label(master, text="Spring constant ( kcal mol"+u"\u207B\u2071"+"%s"%angstrom+u"\u207B\u00B2"+" )" )

# listbox.insert(i, (u"C\u2076"))

        #self.restrain_top2 = Tkinter.Text(master)
        #self.restrain_top2.tag_configure("s", offset=5)
        #self.restrain_top2.insert("insert","Spring constant (kcal mol","","-1","s"," %s"%angstrom,"","-2","s",")")
        #self.restrain_top2.configure(state="disabled")
        #self.restrain_top2.pack(side="top")
        self.restrain_top2.grid(row=14, column=0, sticky=W)
        self.restrain_K = Tkinter.Entry(master, width=9)
        self.restrain_K.insert(0, "2.5")
        self.restrain_K.grid(row=14, column=1, sticky=W)
        Label(master, text="  ").grid(row=15, columnspan=5, sticky=N)

        self.temp_top = Tkinter.Label(master, text="Temperature (K)")
        self.temp_top.grid(row=21, column=0, sticky=W)
        self.temp = Tkinter.Entry(master, width=9)
        self.temp.insert(0, "298.15")
        self.temp.grid(row=21, column=1, sticky=W)

        self.timestep_top = Tkinter.Label(master, text="Timestep (fs)")
        self.timestep_top.grid(row=22, column=0, sticky=W)
        self.timestep = Tkinter.Entry(master, width=9)
        self.timestep.insert(0, "2.0")
        self.timestep.grid(row=22, column=1, sticky=W)

        self.equ_step_top = Tkinter.Label(master, text="Number of Steps for Equilibration")
        self.equ_step_top.grid(row=23, column=0, sticky=W)
        self.equ_step = Tkinter.Entry(master, width=9)
        self.equ_step.insert(0, "1000000")
        self.equ_step.grid(row=23, column=1, sticky=W)

        self.prod_step_top = Tkinter.Label(master, text="Number of Steps for Production")
        self.prod_step_top.grid(row=24, column=0, sticky=W)
        self.prod_step = Tkinter.Entry(master, width=9)
        self.prod_step.insert(0, "10000000")
        self.prod_step.grid(row=24, column=1, sticky=W)

        self.interval_top = Tkinter.Label(master, text="Report Interval")
        self.interval_top.grid(row=25, column=0, sticky=W)
        self.interval = Tkinter.Entry(master, width=9)
        self.interval.insert(0, "1000")
        self.interval.grid(row=25, column=1, sticky=W)

        Label(master, text="  ").grid(row=26, columnspan=5, sticky=N)
        self.dialog.activate(geometry = 'centerscreenalways')

    def searchTopFile(self):
        ftypes=(('prmtop file', '*.prmtop'), ("All files", "*"))
        indir = os.getcwd()
        self.TopFile = askopenfilename(title="Select Amber Topology File", initialdir=indir, filetypes=ftypes)
        if self.TopFile:
            self.ambTop.delete(0, END)
            self.ambTop.insert(0, self.TopFile)
            print self.TopFile

    def searchCrdFile(self):
        ftypes=(('inpcrd file', '*.inpcrd'), ("All files", "*"))
        indir = os.getcwd()
        self.CrdFile = askopenfilename(title="Select Amber Coordinate File", initialdir=indir, filetypes=ftypes)
        if self.CrdFile:
            self.ambCrd.delete(0, END)
            self.ambCrd.insert(0, self.CrdFile)
            print self.CrdFile
        
    def read_base_dir(self):
        indir = os.getcwd()
        base_dir = ""
        base_dir = askdirectory(title="Select prepared project directory", initialdir=indir, mustexist=1)
        self.proj_dir.delete(0, END)
        self.proj_dir.insert(0, base_dir)   

    def writeCmdFile(self):
        curdir = os.getcwd()

        prep_info = open("%s/sys_prep" % proj_dir, "r")
        line = prep_info.readline()
        line = prep_info.readline()
        forcefield = line.strip()
        line = prep_info.readline()
        watermodel = line.strip()
        prep_info.close()

        filename = "%s/run_omm_watsite.sh" % proj_dir
        print filename

        fi = open(filename, "w")
        if forcefield != "AMOEBA-2013":
            fi.write("# prepare the system from amber files\n")
            fi.write("if [ ! -f sys_min.xml ]; then\n")
            fi.write("%s/python %s/bin/openmm_prep.py  -i %s -t %s \\\n" % (python_dir, watsite_dir, ambCrd, ambTop ) )
            if useShake:
                fi.write("-l %s -c %.2f --shake \\\n" % (longrange_method, longrange_cutoff) )
            else:
                fi.write("-l %s -c %.2f \\\n" % (longrange_method, longrange_cutoff) )
            #fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb\n" % restrain_K)
            # make restraint larger during minimization
            if applyRestrain:
                fi.write("--restrain-mask '!:WAT&!@H=' -k 10.0 --reference prot_amber.pdb \\\n")
            fi.write("--cuda %i \n"% (cuda_device))
            fi.write("fi\n\n")

            fi.write("# temperature coupling %.2f K: 25,000 * 1.0 fs = 25 ps\n" % temperature )
            fi.write("if [ ! -f sys_NVT.xml ]; then\n")
            fi.write("%s/python %s/bin/openmm_md.py  -i %s -t %s \\\n" % (python_dir, watsite_dir, ambCrd, ambTop ) )
            fi.write("--xml system.xml -s sys_min.xml --restart sys_NVT.xml \\\n")
            fi.write("-x sys_NVT.nc -r sys_NVT.info -o sys_NVT.out \\\n")
            fi.write("--temp %.2f --gamma_ln 1.0 -n 25000 --interval 1000 --dt 1.0 \\\n" % temperature )
            if applyRestrain:
            	fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb \\\n" % (restrain_K))
            fi.write("--cuda %i \n" % (cuda_device))
            fi.write("fi\n\n")


            fi.write("# pressure coupling %.2f K and 1 bar: %d (stpes) * %.1f (timestep fs) = %.1f ps = %.1f ns\n" % (temperature, equ_step, timestep, equ_step*timestep*0.001, equ_step*timestep*0.000001 ) )
            fi.write("if [ ! -f sys_NPT.xml ]; then\n")
            fi.write("%s/python %s/bin/openmm_md.py  -i %s -t %s \\\n" % (python_dir, watsite_dir, ambCrd, ambTop ) )
            fi.write("--xml system.xml -s sys_NVT.xml --restart sys_NPT.xml \\\n")
            fi.write("-x sys_NPT.nc -r sys_NPT.info -o sys_NPT.out \\\n")
            fi.write("--temp %.2f --gamma_ln 1.0 -n %d --interval 10000 --dt %.1f \\\n" % (temperature, equ_step, timestep ) )
            if applyRestrain:
            	fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb \\\n" % (restrain_K))
            fi.write("--cuda %i --npt\n" % (cuda_device))
            fi.write("fi\n\n")


            fi.write("# production run %.2f K and 1 bar: %d (stpes) * %.1f (timestep fs) = %.1f ps = %.1f ns\n" % (temperature, prod_step, timestep, prod_step*timestep*0.001, prod_step*timestep*0.000001 ) )
            one_md_frames = 500000 # 500,000 steps; x2fs == 1 ns
            split_n = int( prod_step / one_md_frames )
            for i in range(split_n):
                fi.write("if [ ! -f sys_md_%d.xml ]; then\n" %(i+1) )
                fi.write("%s/python %s/bin/openmm_md.py  -i %s -t %s \\\n" % (python_dir, watsite_dir, ambCrd, ambTop ) )
                if i == 0:
                    fi.write("--xml system.xml -s sys_NPT.xml --restart sys_md_%d.xml \\\n" %(i+1) )
                else:
                    fi.write("--xml system.xml -s sys_md_%d.xml --restart sys_md_%d.xml \\\n" %(i, i+1) )
                fi.write("-x sys_md_%d.nc  -r sys_md_%d.info  -o sys_md_%d.out \\\n" % (i+1, i+1, i+1))
                fi.write("--temp %.2f --gamma_ln 1.0 -n %d --interval %d --dt %.1f \\\n" % (temperature, one_md_frames, interval, timestep ) )
                
            	if applyRestrain:
            	    fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb \\\n" % (restrain_K))
                fi.write("--water %s --prev_frame %d -p prot_amber.pdb --cuda %i --npt \n" %(watermodel, one_md_frames/interval*i, cuda_device) )
                fi.write("fi\n\n")
            fi.write("\n")
        else:
            fi.write("# create the system from PDB file, and assign AMOEBA force field parameter in OpenMM\n")
            if watermodel == "AMOEBA03":
                fi.write("%s/python %s/bin/openmm_prep.py  -p prot_amber.pdb -f amoeba2013.xml \\\n" % (python_dir, watsite_dir) )
            #elif watermodel == "AMOEBA14":
            #    fi.write("%s/python %s/bin/openmm_prep.py  -p prot_amber.pdb -f iamoeba.xml amoeba2013xml \\\n" % (python_dir,watsite_dir) )
            else:
            #    tkMessageBox.showwarning("Warning","Choose from two water models AMOEBA03 or AMOEBA14 for AMOEBA-2013 force field...\n")
                tkMessageBox.showwarning("Warning","Only AMOEBA03 water model is compatible with AMOEBA-2013 force field...\n")
                #self.dialog.deactivate()
                self.dialog.withdraw()
                return -1
            fi.write("-l %s -c %.2f -v %.2f --cuda %i  \\\n" % (longrange_method, longrange_cutoff, (longrange_cutoff+2.0), cuda_device) )
            fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb\n" % restrain_K)
            fi.write("\n")

            fi.write("# temperature coupling %.2f K: 25,000 * 1.0 fs = 25 ps\n" % temperature )
            fi.write("%s/python %s/bin/openmm_md.py  -p prot_amber.pdb \\\n" % (python_dir, watsite_dir) )
            fi.write("--xml system.xml -s sys_min.xml --restart sys_NVT.xml \\\n")
            fi.write("-x sys_NVT.nc -r sys_NVT.info -o sys_NVT.out \\\n")
            fi.write("--temp %.2f --gamma_ln 1.0 -n 25000 --interval 2500 --dt 1.0 \\\n" % temperature )
            fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb --cuda %i\n" % (restrain_K, cuda_device) )
            fi.write("\n")

            fi.write("# pressure coupling %.2f K and 1 bar: %d (stpes) * %.1f (timestep fs) = %.1f ps = %.1f ns\n" % (temperature, equ_step, timestep, equ_step*timestep*0.001, equ_step*timestep*0.000001 ) )
            fi.write("%s/python %s/bin/openmm_md.py  -p prot_amber.pdb \\\n" % (python_dir,watsite_dir) )
            fi.write("--xml system.xml -s sys_NVT.xml --restart sys_NPT.xml \\\n")
            fi.write("-x sys_NPT.nc -r sys_NPT.info -o sys_NPT.out \\\n")
            fi.write("--temp %.2f --gamma_ln 1.0 -n %d --interval 10000 --dt %.1f \\\n" % (temperature, equ_step, timestep ) )
            fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb --cuda %i --npt\n" % (restrain_K, cuda_device) )
            fi.write("\n")

            fi.write("# production run %.2f K and 1 bar: %d (stpes) * %.1f (timestep fs) = %.1f ps = %.1f ns\n" % (temperature, prod_step, timestep, prod_step*timestep*0.001, prod_step*timestep*0.000001 ) )
            fi.write("%s/python %s/bin/openmm_md.py  -p prot_amber.pdb \\\n" % (python_dir, watsite_dir) )
            fi.write("--xml system.xml -s sys_NPT.xml --restart sys_fin.xml \\\n")
            fi.write("-x sys_md.nc  -r sys_md.info  -o sys_md.out \\\n")
            fi.write("--temp %.2f --gamma_ln 1.0 -n %d --interval %d --dt %.1f \\\n" % (temperature, prod_step, interval, timestep ) )
            fi.write("--restrain-mask '!:WAT&!@H=' -k %.1f --reference prot_amber.pdb \\\n" % restrain_K)
            fi.write("--water AMOEBA --cuda %i --npt \n" %(cuda_device) )
            fi.write("\n")
            fi.write("\n")

        #"--restrain '!:WAT&!@H=' -k 5 --reference prot_amber.pdb"

        os.system("chmod 777 %s" % filename)

        fmon = open("%s/WATsite.out" % (proj_dir), "w")
        fmon.write("ProductionStep      %i\n" % (prod_step) )
        fmon.write("Interval            %i\n" % (interval) )
        fmon.write("numFrame            %i\n" % (prod_step / interval) )
        fmon.write("Protein             MyProtein.pdb\n")
        if os.path.isfile("MyLigand.mol2"):
            fmon.write("Ligand              MyLigand.mol2\n")
        fmon.write("HydrationSiteMol2   WATsite_OUT/HydrationSites.mol2\n")
        fmon.write("HydrationSitePDB    WATsite_OUT/HydrationSites.pdb\n")
        fmon.write("GridMol2            WATsite_OUT/FilterGrid.mol2\n")
        fmon.write("GridDX              WATsite_OUT/grid_occupancy.dx\n")
        fmon.write("WaterTraj           WATsite_OUT/WATinside.mol2\n")
        fmon.write("EnergyFile          WATsite_OUT/cluster.egy\n")
        fmon.close()

    def apply(self, result):
        global proj_dir
        global ambTop
        global ambCrd
        global cuda_device
        global longrange_method
        global longrange_cutoff
        global useShake
        global applyRestrain
        global restrain_K
        global temperature
        global timestep
        global equ_step
        global prod_step
        global interval
        
        if result == 'Cancel':
            self.dialog.deactivate()
            self.dialog.withdraw()
        elif result == 'Save Parameters':
            proj_dir         = self.proj_dir.get()
            ambTop           = self.ambTop.get()
            ambCrd           = self.ambCrd.get()
            cuda_device      = int(self.cuda_index.get())
            longrange_method = self.lr.get()
            longrange_cutoff = float(self.lr_cut.get())
            useShake         = int(self.shake.get())
            applyRestrain    = int(self.restrain.get())
            restrain_K       = float(self.restrain_K.get())
            temperature      = float(self.temp.get())
            timestep         = float(self.timestep.get())
            equ_step         = int(self.equ_step.get())
            prod_step        = int(self.prod_step.get())
            interval         = int(self.interval.get())

            if not os.path.isfile("prot_amber.pdb"):
                command = '%s/bin/ambpdb -p %s -c %s -ep > prot_amber.pdb' % (amber_home, ambTop, ambCrd)
                print "Executing ambpdb...\n %s" % command
                os.system(command)

            self.writeCmdFile()
            tkMessageBox.showinfo("Preparation Done","1. Enter folder: %s\n2. Start OpenMM simulation by 'nohup ./run_openmm.sh &'.\n" % (self.proj_dir.get()))

            self.dialog.deactivate()
            self.dialog.withdraw()

def setParam4MD(app):
    setSimulationParameter(app.root)


#=====================================================================================#
#=======================    Menu 4: Perform WATsite Analysis    ======================#
#=====================================================================================#
class performWATsiteAnalysis:

    def __init__(self, top):

        fi = open("sys_prep","r")
        line = fi.readline()
        line = fi.readline()
        read_watermodel = fi.readline().strip()
        fi.close()

        fin = open("WATsite.out", "r")
        for line in fin:
            if line.find("numFrame") > -1:
                numFrame = float(line.split()[1])
            if line.find("ProductionStep") > -1:
                prod_step = float(line.split()[1])
            if line.find("Interval") > -1:
                interval = float(line.split()[1])
        fin.close()

        # Create the dialog.
        self.dialog = Pmw.Dialog(top,
                                 buttons = ('Perform Analysis', 'Cancel'),
                                 defaultbutton = 'Cancel',
                                 title = 'PyMOL WATsite Plugin',
                                 command = self.apply)
        self.curdir = os.getcwd()
        
        w1 = Tkinter.Label(self.dialog.interior(),
                          text = MAL_INFO, 
                          background = 'black',
                          foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='we')
        w1.grid_columnconfigure(0, weight=1)

        master = self.dialog.interior()

        Label(master, text='\nPerform WATsite Analysis\n', foreground = 'blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=5)

        curdir = os.getcwd()
        Tkinter.Label(master, text="Prepared Amber Topology File").grid(row=3, column=0, sticky=W)
        self.ambTop = Tkinter.Entry(master, width=50)
        self.ambTop.insert(0, "%s/prot_amber.prmtop" % curdir)
        self.ambTop.grid(row=3, column=1)
        self.butTop = Tkinter.Button(master, text='Browse', command=self.searchTopFile )
        self.butTop.grid(row=3, column=2)

        Tkinter.Label(master, text="Prepared Amber Coordinate File").grid(row=4, column=0, sticky=W)
        self.ambCrd = Tkinter.Entry(master, width=50)
        self.ambCrd.insert(0, "%s/prot_amber.inpcrd" % curdir)
        self.ambCrd.grid(row=4, column=1)
        self.butCrd = Tkinter.Button(master, text='Browse', command=self.searchCrdFile )
        self.butCrd.grid(row=4, column=2)

        Tkinter.Label(master, text="Simulation Trajectory File").grid(row=5, column=0, sticky=W)
        self.TrajFile = Tkinter.Entry(master, width=50)
        self.TrajFile.insert(0, "%s/sys_md_?.nc" % curdir)
        self.TrajFile.grid(row=5, column=1)
        self.butTraj = Tkinter.Button(master, text='Browse', command=self.searchNcFile )
        self.butTraj.grid(row=5, column=2)

        Tkinter.Label(master, text="  ").grid(row=10, columnspan=5, sticky=N)

        """
        self.timestep_top = Tkinter.Label(master, text="Timestep (fs)")
        self.timestep_top.grid(row=11, column=0, sticky=W)
        self.timestep = Tkinter.Entry(master, width=9)
        self.timestep.insert(0, "1.0")
        self.timestep.grid(row=11, column=1, sticky=W)
        """

        self.prod_step_top = Tkinter.Label(master, text="Number of Steps for Production")
        self.prod_step_top.grid(row=12, column=0, sticky=W)
        self.prod_step = Tkinter.Entry(master, width=9)
        self.prod_step.insert(0, "%d" % prod_step)
        self.prod_step.grid(row=12, column=1, sticky=W)

        self.interval_top = Tkinter.Label(master, text="Report Interval")
        self.interval_top.grid(row=13, column=0, sticky=W)
        self.interval = Tkinter.Entry(master, width=9)
        self.interval.insert(0, "%d" % interval)
        self.interval.grid(row=13, column=1, sticky=W)

        Tkinter.Label(master, text="  ").grid(row=20, columnspan=5, sticky=N)

        self.wat_top1 = Tkinter.Label(master, text="Water Model")
        self.wat_top1.grid(row=21, column=0, sticky=W)
        self.wat = Tkinter.StringVar()
        self.wat.set(read_watermodel)
        self.wat_top2 = Tkinter.OptionMenu(master, self.wat, 'SPC/E', 'TIP3P', 'OPC', 'TIP4P', 'TIP4PEW', 'AMOEBA03')
        self.wat_top2.grid(row=21, column=1, sticky=W)
      
        self.cluster_top1 = Tkinter.Label(master, text="Clustering Method")
        self.cluster_top1.grid(row=22, column=0, sticky=W)
        self.cluster = Tkinter.StringVar()
        self.cluster.set("QT")
        self.cluster_top2 = Tkinter.OptionMenu(master, self.cluster, 'DBSCAN', 'QT')
        self.cluster_top2.grid(row=22, column=1, sticky=W)

        Label(master, text="  ").grid(row=30, columnspan=5, sticky=N)
        self.dialog.activate(geometry = 'centerscreenalways')

    def searchTopFile(self):
        ftypes=(('prmtop file', '*.prmtop'), ("All files", "*"))
        indir = os.getcwd()
        self.TopFile = askopenfilename(title="Select Amber Topology File", initialdir=indir, filetypes=ftypes)
        if self.TopFile:
            self.ambTop.delete(0, END)
            self.ambTop.insert(0, self.TopFile)
    def searchCrdFile(self):
        ftypes=(('inpcrd file', '*.inpcrd'), ("All files", "*"))
        indir = os.getcwd()
        self.CrdFile = askopenfilename(title="Select Amber Coordinate File", initialdir=indir, filetypes=ftypes)
        if self.CrdFile:
            self.ambCrd.delete(0, END)
            self.ambCrd.insert(0, self.CrdFile)
    def searchNcFile(self):
        ftypes=(('netcdf file', '*.nc'), ("All files", "*"))
        indir = os.getcwd()
        self.Traj = askopenfilename(title="Select Trajectory File", initialdir=indir, filetypes=ftypes)
        if self.Traj:
            self.TrajFile.delete(0, END)
            self.TrajFile.insert(0, self.Traj)
            print self.Traj

    def align_amber_crystal(self):

        if not os.path.isfile("prot_amber.pdb"):
            command = '%s/bin/ambpdb -p prot_amber.prmtop -c prot_amber.inpcrd -ep > prot_amber.pdb' % (amber_home)
            print "Executing ambpdb...\n %s" % command
            os.system(command)

        cmd.delete("all")
        cmd.load("prot_amber.pdb")
        cmd.load("MyProtein.pdb")
        cmd.super("prot_amber", "MyProtein")
        cmd.save("ref.pdb", "prot_amber", 0, "pdb")

    def runCpptraj(self):
        curdir = os.getcwd()
        filename = "%s/cpptraj.in" % curdir

        fi = open(filename, "w")
        if watermodel == 'AMOEBA03':
            fi.write("parm %s/prot_amber.pdb\n" % curdir)
        else:
            fi.write("parm %s [sys]\n" % ambTop)
        fi.write("parm %s/ref.pdb [ref]\n" % curdir)
        fi.write("reference %s/ref.pdb parm [ref]\n" % curdir)
        fi.write("\n")

        one_md_frames = 500000 # 500,000 steps; x2fs == 1 ns
        split_n = int( prod_step / one_md_frames )

        if watermodel == 'AMOEBA03':
            fi.write("trajin %s parm %s\n" % (TrajFile, 'prot_amber.pdb'))
        else:
            for i in range(1, split_n+1):
                fi.write("trajin sys_md_%d.nc parm [sys]\n" % (i) )
        fi.write("\nautoimage\n")
        fi.write("rms reference !(:WAT,Cl\-,Na\+,CA)&!@H= out rms.dat\n")
        fi.write("\n")
        fi.write("trajout mdsnaps.pdb pdb include_ep\n")
        fi.write("go\n")
        fi.write("\n")
        fi.write("rms reference !(:WAT,Cl\-,Na\+,CA)&!@H= out rms.dat\n")
        fi.write("strip :Cl-\n")
        fi.write("strip :Na+\n")
        fi.write("strip :WAT,HOH\n")
        fi.write("trajout firstFrame.pdb pdb onlyframes 1 sp include_ep\n")
        fi.write("\n")
        fi.write("go\nquit\n")
        fi.write("\n")
        fi.close()
        print "\n\n"+"#"+"="*25+"#"
        command = '%s/bin/cpptraj < cpptraj.in' % (amber_home)
        print "Executing cpptraj...\n %s" % command
        os.system(command)
        print "#"+"="*25+"#\n\n"

    def writeBcfFile(self):
        curdir = os.getcwd()
        filename = "%s/input.bcf" % curdir
        print filename

        if not os.path.isfile("MyBindingSite.pdb"):
            tkMessageBox.showwarning(
            "Check Input",
            "Cannot find binding site file: MyBindingSite.pdb")
            return -1

        if not os.path.isfile("mdsnaps.pdb"):
            tkMessageBox.showwarning(
            "Check Input",
            "Cannot find trajectory file: mdsnaps.pdb")
            return -1

        if not os.path.isdir("waterEnergy"):
            tkMessageBox.showwarning(
            "Check Input",
            "Cannot find water energy directory: waterEnergy")
            return -1

        if not 'margin' in globals():
            fi = open("sys_prep","r")
            margin = float(fi.readline().strip())
            fi.close()
        if not margin:
            margin = 2.5 #default
        print "Create 3D grid box by extending %.2f A from binding site structure" % margin

        fi = open(filename, "w")
        fi.write("$Input_files:\n")
        fi.write("%20s          %-50s\n" % ("MyBindingSite.pdb", "|input pdb file to define protein binding site")  )
        fi.write("%20s          %-50s\n" % ("../mdsnaps.pdb", "|input trajectory in pdb")  )
        fi.write("%20s          %-50s\n" % ("../waterEnergy", "|folder contain water energy files")  )
        fi.write("\n")

        fi.write("$Input_MDsnapshots:\n")
        fi.write("%20i          %-50s\n" % (int(numFrame), "|number of frames in trajectory files")  )
        if watermodel == 'AMOEBA03':
            fi.write("%20s          %-50s\n" % ('AMOEBA', "|water model (SPC/E; TIP3P; OPC; TIP4P; TIP4PEW; AMOEBA)")  )
        else:
            fi.write("%20s          %-50s\n" % (watermodel, "|water model (SPC/E; TIP3P; OPC; TIP4P; TIP4PEW; AMOEBA)")  )
        fi.write("\n")

        fi.write("$Grid_Parameters:\n")
        fi.write("%20s          %-50s\n" % ("1.0", "|size of water (SD)")  )
        fi.write("%20s          %-50s\n" % ("0.25", "|griddelta: distance between two adjacent grid points")  )
        fi.write("%20s          %-50s\n" % (margin, "|Farest Distance to extend binding site box")  )
        fi.write("%20s          %-50s\n" % ("0.045", "|Water density for grid cutoff")  )
        fi.write("\n")

        fi.write("$Cluster_method:\n")
        if cluster_method == "DBSCAN":
            fi.write("%20s          %-50s\n" % ("1", "|clustering method for HS identification (=1: DBSCAN; =2: QT clustering)")  )
        else:
            fi.write("%20s          %-50s\n" % ("2", "|clustering method for HS identification (=1: DBSCAN; =2: QT clustering)")  )
        fi.write("\n")

        fi.write("$DBSCAN:\n")
        fi.write("%20s          %-50s\n" % ("300", "|min of points define neighbor in DBSCAN (start number)")  )
        fi.write("%20s          %-50s\n" % ("80", "|min of points define neighbor in DBSCAN (end number)")  )
        fi.write("\n")

        fi.write("$Hydro_Cluster:\n")
        fi.write("%20s          %-50s\n" % ("10", "|maximum clusters")  )
        fi.write("%20s          %-50s\n" % ("2.50", "|cluster_mean_dist")  )
        fi.write("%20s          %-50s\n" % ("2.75", "|maxdist")  )
        fi.write("%20s          %-50s\n" % ("10.0", "|distance cutoff for QT clustering")  )
        fi.write("\n")

        fi.write("$Entropy:\n")
        fi.write("%20s          %-50s\n" % ("9", "|covariance dimension (=3: two 3x3 matrix; =6: one 6x6 dimension; =9: both 3x3 and 6x6 method)")  )
        fi.write("%20s          %-50s\n" % ("70", "|number of bins used to construct Probability Distribution function")  )
        fi.write("\n")
        return 0

    def apply(self, result):

        global ambTop
        global ambCrd
        global TrajFile
        #global timestep
        global prod_step
        global interval
        global watermodel
        global cluster_method
        global numFrame

        if result == 'Cancel':
            self.dialog.deactivate()
            self.dialog.withdraw()
        elif result == 'Perform Analysis':
            ambTop           = self.ambTop.get()
            ambCrd           = self.ambCrd.get()
            TrajFile         = self.TrajFile.get()
            #timestep         = float(self.timestep.get())
            prod_step        = int(self.prod_step.get())
            interval         = int(self.interval.get())
            watermodel       = self.wat.get()
            cluster_method   = self.cluster.get()

            if not os.path.isfile("prot_amber.pdb"):
                command = '%s/bin/ambpdb -p %s -c %s -ep > prot_amber.pdb' % (amber_home, ambTop, ambCrd)
                print "Executing ambpdb...\n %s" % command
                os.system(command)

            numFrame = int(prod_step*1.0/interval)
            self.align_amber_crystal()
            if not os.path.isfile("mdsnaps.pdb"):
                self.runCpptraj()

            flag = self.writeBcfFile()
            if flag == 0:
                flag_wat = os.system("%s/bin/watsite -c input.bcf -o WATsite_OUT" % watsite_dir)
                if flag_wat == 0:
                    tkMessageBox.showinfo("WATsite analysis done","Investigate results by 'Import WATsite Results'.")

                    self.dialog.deactivate()
                    self.dialog.withdraw()

def runWATsite(app):
    performWATsiteAnalysis(app.root)



#=======================================================================================#
#=====================   Menu 5: Import & Display WATsite Results   ====================#
#=======================================================================================#
def import_results(app):

    global import_dir
    global HydrationSiteEnergyFile
    global HydrationSiteMol2
    global prot_obj
    global lig_obj
    global HS_mol2_obj
    global num_WATsite
    global TdS
    global dH
    global dG
    global Occupancy
    global numFrame
    numFrame = 5000.0

    global grid_mol2, grid_dx

    ftypes=(('out file', '*.out'), ("All files", "*"))
    curdir = os.getcwd()
    openfile = askopenfilename(initialdir=curdir, filetypes=ftypes)

    #print "\n\n\n"
    #print openfile

    if openfile:
        import_dir = os.path.dirname(openfile)
        os.chdir(import_dir)
        
        fin = open(openfile, "r")
        HSwithLig = 0
        for line in fin:
            if line.find("Protein") > -1:
                ProtFile = line.split()[1]
                prot_obj = ProtFile.split("/")[-1].split(".")[0]
            if line.find("Ligand") > -1:
                LigFile = line.split()[1]
                lig_obj = LigFile.split("/")[-1].split(".")[0]
                HSwithLig = 1
            if line.find("HydrationSiteMol2") > -1:
                HydrationSiteMol2 = line.split()[1]
                HS_mol2_obj = HydrationSiteMol2.split("/")[-1].split(".")[0]
            if line.find("HydrationSitePDB") > -1:
                HydrationSitePDB = line.split()[1]
            if line.find("GridMol2") > -1:
                Grid_mol2 = line.split()[1]
                Grid_mol2_obj = 'Grid_mol2'
            if line.find("GridDX") > -1:
                Grid_dx = line.split()[1]
                Grid_dx_obj = 'Grid_dx'
            if line.find("WaterTraj") > -1:
                WaterTraj = line.split()[1]
            if line.find("EnergyFile") > -1:
                HydrationSiteEnergyFile = line.split()[1]
            if line.find("numFrame") > -1:
                numFrame = float(line.split()[1])
        fin.close()

        #print "Energy File is: %s" % HydrationSiteEnergyFile 

        flag = 0
        if os.path.isfile(HydrationSiteEnergyFile):
            flag = 1
        if flag == 0:
            tkMessageBox.showwarning(
                "Could not find Hydratioin Site Energy files",
                "WATsite might be still in progress or didn't finish successfully"
            )
            return 0

        TdS, dH, dG, Occupancy = [], [], [], []
        sie_file = open(HydrationSiteEnergyFile, 'r')
        for i in sie_file:
            if i.find("@") == -1 and i != "":
                j = i.split()
                S_val = float(j[1])
                H_val = float(j[2])
                TdS.append(S_val)
                dH.append(H_val)
                dG.append(S_val+H_val)
                Occupancy.append(float(j[8]))
        sie_file.close()
        num_WATsite = len(dG)
        z = loadDialog(app.root)

        if z.ok_flag == 0:
            return
        if tkMessageBox.askokcancel("Read project", "Reading results will delete all exisiting objects in current session.") == 0:
            return 0
        # remove existing objects
        for na in cmd.get_names("objects"):
            try:
                cmd.remove(na)
            except:
                pass
            cmd.delete(na)

        if z.protfile.get() == 1:
            try:
                fi = open(ProtFile, 'r')
            except:
                tkMessageBox.showwarning(
                "Missing protein file",
                "Please, check the protein file specified in WATsite.out is existing in the directory"
                )
                return 0
            cmd.load(ProtFile)
            cmd.show("cartoon", prot_obj)
            util.cbc(prot_obj)
            util.cnc(prot_obj)
        if z.ligfile.get() == 1 and HSwithLig:
            try:
                fi = open(LigFile, 'r')
            except:
                tkMessageBox.showwarning(
                "Missing ligand file",
                "Please, check the ligand file specified in WATsite.out is existing in the directory"
                )
                return 0
            cmd.load(LigFile)
            cmd.show("stick", lig_obj)
            cmd.set("stick_radius", "0.15")
            cmd.remove("%s and hydrogen"%lig_obj)
            util.cbaw(lig_obj)

        if z.HSmol2file.get() == 1:
            try:
                fi = open(HydrationSiteMol2, 'r')
            except:
                tkMessageBox.showwarning(
                "Missing hydration site file",
                "Please, check the hydration site mol2 file specified in WATsite.out is existing in the directory"
                )
                return 0

            cmd.load(HydrationSiteMol2)
            cmd.show_as("spheres",HS_mol2_obj)
            cmd.set("sphere_scale", 0.35, HS_mol2_obj)
            cmd.spectrum("pc", "blue_white_red", HS_mol2_obj, "-3","3")
            cmd.label("%s"%HS_mol2_obj,"ID")
            cmd.set("label_size", -1)
            if HSwithLig:
                cmd.distance("hb","%s or %s or (%s and not resn HOH) " % (HS_mol2_obj, lig_obj, prot_obj),"%s or %s or (%s and not resn HOH) " % (HS_mol2_obj, lig_obj, prot_obj), 3.2, quiet=1,mode=2,label=0,reset=1)
            else:
                cmd.distance("hb","%s or (%s and not resn HOH) " % (HS_mol2_obj, prot_obj),"%s or (%s and not resn HOH) " % (HS_mol2_obj, prot_obj), 3.2, quiet=1,mode=2,label=0,reset=1)

        if z.DensityMol2File.get() == 1:
            try:
                fi = open(Grid_mol2, 'r')
            except:
                tkMessageBox.showwarning(
                "Missing grid mol2 file",
                "Please, check the GridMol2 file specified in WATsite.out is existing in the directory"
                )
                return 0
            cmd.load(Grid_mol2, Grid_mol2_obj)
            cmd.set("sphere_scale",0.025,selection="%s"%Grid_mol2_obj)
            cmd.show_as("spheres",Grid_mol2_obj)
            #cmd.hide("nonbonded", "tmp")
            cmd.spectrum("pc", "blue_white_red", Grid_mol2_obj, "0.045","0.2")

        if z.DensityDXFile.get() == 1:
            try:
                fi = open(Grid_dx, 'r')
            except:
                tkMessageBox.showwarning(
                "Missing grid dx file",
                "Please, check the GridDX file specified in WATsite.out is existing in the directory"
                )
                return 0
            cmd.load(Grid_dx, Grid_dx_obj)
            cmd.isomesh("WaterDensity", Grid_dx_obj, level=0.045)
            cmd.color("gray50", "WaterDensity")
            cmd.set("mesh_width", 0.5)
        cmd.set("ray_trace_fog", 0)
        cmd.set("ray_shadows", 0)

        HydrationSiteEnergy(app.root)

class loadDialog(tkSimpleDialog.Dialog):
    def body(self, master):
        self.master.title('PyMOL WATsite Plugin')
        self.var_check = []

        self.ok_flag = 0

        Label(master, text="Please select files to load").grid(row=0, column=0, sticky=W)
        self.protfile = IntVar()
        self.a = Checkbutton(master, text="Protein", variable=self.protfile)
        self.a.grid(row=1, column=0, sticky=W, padx=5)
        self.protfile.set(1)
        
        self.ligfile = IntVar()
        self.b = Checkbutton(master, text="Ligand", variable=self.ligfile)
        self.b.grid(row=2, column=0, sticky=W, padx=5)
        self.ligfile.set(0)
        
        self.HSmol2file = IntVar()
        self.c = Checkbutton(master, text="Hydration Site", variable=self.HSmol2file)
        self.c.grid(row=3, column=0, sticky=W, padx=5)
        self.HSmol2file.set(1)
        
        self.DensityMol2File = IntVar()
        self.c = Checkbutton(master, text="Density Mol2", variable=self.DensityMol2File)
        self.c.grid(row=4, column=0, sticky=W, padx=5)
        self.DensityMol2File.set(1)
        
        self.DensityDXFile = IntVar()
        self.c = Checkbutton(master, text="Density DX", variable=self.DensityDXFile)
        self.c.grid(row=5, column=0, sticky=W, padx=5)
        self.DensityDXFile.set(1)
        
    def apply(self):
        self.ok_flag = 1

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

class HydrationSiteEnergy:
    def __init__(self, top):
        self.dialog = Pmw.Dialog(top, 
                                 buttons = ('OK','Cancel'), 
                                 defaultbutton = 'OK', 
                                 title = 'WATsite Results', 
                                 command = self.apply)        
        parent = self.dialog.interior()

        group1 = Pmw.Group(parent, tag_text='Energy values in kcal/mol (double click to select hydration site in Pymol)')
        master = Frame(group1.interior())
        scrollbar = Scrollbar(master)
        listbox = Tkinter.Listbox(master, yscrollcommand=scrollbar.set)
        listbox.config(font=('Courier', 12))

        TXT='HS#      -T%sS      %sH        %sG    occupancy'%(delta, delta, delta)
        listbox.insert(END, TXT)
        sie_file = open(HydrationSiteEnergyFile, 'r')
        ci = 3
        hsi = 1
        pos = 0

        self.S, self.H, self.G, self.O = [], [], [], []
        for i in sie_file:
            if i.find("@") == -1 and i != "":
                j = i.split()
                S_val = float(j[1])
                H_val = float(j[2])
                G_val = S_val+H_val
                O_val = float(j[8])
                self.S.append(S_val)
                self.H.append(H_val)
                self.G.append(G_val)
                self.O.append(O_val)
                TXT = '%3d    %6.2f    %6.2f    %6.2f    %6.2f' % (hsi, S_val, H_val, G_val, O_val)
                listbox.insert(END, TXT)
                #pos += 1
                ci += 1
                hsi += 1
        sie_file.close()
        
        scrollbar.config(command=listbox.yview)
        listbox.bind('<Double-1>', self.handleList)
        self.listbox = listbox
        scrollbar.pack(side=RIGHT, fill=Y)
        listbox.pack(side=LEFT, expand=YES, fill=BOTH)
        master.pack(fill = 'both', expand = 1, padx = 5, pady=1)
        group1.pack(fill = 'both', expand = 1, padx = 5, pady = 5)
        
        
        group2 = Pmw.Group(parent, tag_text='Select how to color the hydration sites')
        f2 = Frame(group2.interior())
        
        b1 = Button(f2, text='-T%sS'%delta, command=self.colorByS)
        b1.grid(row=2, column=2)
        b2 = Button(f2, text='%sH'%delta, command=self.colorByH)
        b2.grid(row=2, column=4)
        b3 = Button(f2, text='%sG'%delta, command=self.colorByG)
        b3.grid(row=2, column=5)
        b4 = Button(f2, text='Occupancy', command=self.colorByO)
        b4.grid(row=2, column=6)
        f2.pack(fill = 'both', expand = 0, padx = 5, pady=1)
        group2.pack(fill = 'both', expand = 0, padx = 5, pady = 5)
        
        #self.dialog.activate(geometry = 'centerscreenalways')

    def colorByS(self):
        cmd.select("colorByS", HS_mol2_obj)
        for i in range(len(self.S)):
            cmd.alter("colorByS and id %s"%(i+1), "b=%s"%self.S[i])
        #cmd.spectrum("b", "blue_white_red", "colorByS",minimum=min(self.S), maximum=max(self.S))
        cmd.spectrum("b", "blue_white_red", "colorByS",minimum=0, maximum=2)
        cmd.delete("colorByS")
    def colorByH(self):
        cmd.select("colorByH", HS_mol2_obj)
        for i in range(len(self.H)):
            cmd.alter("colorByH and id %s"%(i+1), "b=%s"%self.H[i])
        #cmd.spectrum("b", "blue_white_red", "colorByH",minimum=min(self.H), maximum=max(self.H))
        cmd.spectrum("b", "blue_white_red", "colorByH",minimum=-5, maximum=5)
        cmd.delete("colorByH")
    def colorByG(self):
        cmd.select("colorByG", HS_mol2_obj)
        for i in range(len(self.G)):
            cmd.alter("colorByG and id %s"%(i+1), "b=%s"%self.G[i])
        #cmd.spectrum("b", "blue_white_red", "colorByG",minimum=min(self.G), maximum=max(self.G))
        cmd.spectrum("b", "blue_white_red", "colorByG",minimum=-3, maximum=3)
        cmd.delete("colorByG")
    def colorByO(self):
        cmd.select("colorByO", HS_mol2_obj)
        for i in range(len(self.O)):
            cmd.alter("colorByO and ID %s"%(i+1), "b=%s"%self.O[i])
        #cmd.spectrum("b", "blue_white_red", "colorByO",minimum=min(self.G), maximum=max(self.G))
        cmd.spectrum("b", "blue_white_red", "colorByO",minimum=0.3, maximum=1.0)
        cmd.delete("colorByO")           
    def handleList(self, event):
        index = self.listbox.curselection()
        label = self.listbox.get(index)
        self.runCommand(label)
    def runCommand(self, selection):
        hs_index = selection.split()[0]
        cmd.select("selected_hs", HS_mol2_obj + ' and id '+hs_index)
        cmd.zoom("selected_hs")
        #print HS_mol2_obj, hs_index
    def apply(self, result):
        if result == 'OK':
            self.ok_flag = 1
            self.dialog.deactivate()
            self.dialog.withdraw()
        else:
            self.ok_flag = 0
            self.dialog.deactivate()
            self.dialog.withdraw()




#=======================================================================================#
#=================== Menu 6: Estimate Protein Desolvation Free Energy ==================#
#=======================================================================================#
def estimate_desolvation(app):       
    global liglist
    global num_ligands
    global desol_G, desol_S, desol_H #desolvation energies
    desol_G, desol_S, desol_H = [], [], []
    
    imp_lig = import_ligands(app.root)
    
    sel_dir = imp_lig.project_dir
    liglist = []
    for lf in os.listdir(sel_dir):
        if lf.find("mol2") > -1 or lf.find("pdb") > -1:
            cmd.load("%s/%s"%(sel_dir,lf))
            ligname = lf.split(".")[0]
            cmd.remove("%s and hydrogen"%ligname)
            cmd.show_as("sticks", ligname)
            liglist.append(ligname)
    num_ligands = len(liglist)
    cmd.set("stick_radius", "0.15")
    try:
        HydrationSiteMol2
    except NameError:
        tkMessageBox.showwarning(
        "Bad input",
        "Please import hydration site results first"
        )
        return 0

    hs = HS_mol2_obj

    if imp_lig.useScale:
        disCutoff = imp_lig.HSSradius
        for lig in liglist:
            sel_hs = '%s_replacedHS' % lig
            cmd.select(sel_hs, '%s w. %.2f of %s'%(hs, disCutoff, lig))

            ligG, ligS, ligH = 0.0, 0.0, 0.0
            alldis = []
            numLigAtoms = cmd.count_atoms('%s' % lig)

            for j in range(num_WATsite):
                for i in range(numLigAtoms):
                    dis = cmd.get_distance(atom1='%s and idx %d' % (lig, i+1), atom2='%s and id %d' % (hs, j+1) )
                    alldis.append(dis)
                min_dis = min(alldis)
                if min_dis < disCutoff:
                    ligG += (0 - dG[j] )*(1-min_dis/disCutoff)
                    ligS += (0 - TdS[j])*(1-min_dis/disCutoff)
                    ligH += (0 - dH[j] )*(1-min_dis/disCutoff)
            print lig, ligG, ligS, ligH
            desol_G.append(ligG)
            desol_S.append(ligS)
            desol_H.append(ligH)

    else:
        #set b factor values by TdS, dH, dG
        for i in range(num_WATsite):
            cmd.alter('%s and id %s'%(hs,i+1), "b=%s"%dG[i])
            cmd.alter('%s and id %s'%(hs,i+1), "partial_charge=%s"%TdS[i])
            cmd.alter('%s and id %s'%(hs,i+1), 'vdw=%s'%dH[i])
        ligG, ligS, ligH = [], [], []
        for lig in liglist:
            sel_hs = '%s_replacedHS'%lig
            cmd.select(sel_hs, '%s w. %.2f of %s'%(hs, imp_lig.HSSradius, lig))
            #print '%s w. %.2f of %s'%(hs, imp_lig.HSSradius, lig)
            myspace = {'ligG': [], 'ligS': [], 'ligH': []}
            cmd.iterate(sel_hs, 'ligG.append(b)', space=myspace)
            cmd.iterate(sel_hs, 'ligS.append(partial_charge)', space=myspace)
            cmd.iterate(sel_hs, 'ligH.append(vdw)', space=myspace)
            #print myspace['ligG']
            sumG, sumS, sumH = 0.0, 0.0, 0.0
            for i in myspace['ligG']:
                    sumG += (0 - float(i))
            for i in myspace['ligS']:
                    sumS += (0 - float(i))
            for i in myspace['ligH']:
                    sumH += (0 - float(i))
            desol_G.append(sumG)
            desol_S.append(sumS)
            desol_H.append(sumH)
            print lig, sumG, sumS, sumH
    displayResultsDialog(app.root)
        
class import_ligands(tkSimpleDialog.Dialog):
    def body(self, master):
        self.master.title('PyMOL WATsite Plugin')
        self.applic = master
        self.ok_flag = 0

        w1 = Tkinter.Label(master,
                          text = MAL_INFO,
                          background = 'black',
                          foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='we')
        w1.grid_columnconfigure(0, weight=1)

        Label(master, text='\nEstimate Desolvation Energy for Ligand(s)\n', foreground = 'blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=5)

        curdir = os.getcwd()
        Label(master, text="Directory Contains All Ligands").grid(row=2, column=1, sticky=W)

        self.eproj = Entry(master, width=35)
        self.eproj.insert(0, "Full path to ligand library")
        self.eproj.grid(row=2, column=2)
        self.bproj = Button(master, text='Browse', command=self.read_base_dir)
        self.bproj.grid(row=2, column=3)
        
        Label(master, text="Hydration Site Selection Radius (%s)" % angstrom).grid(row=3, column=1, sticky=W)

        self.hssr = Entry(master, width=10)
        self.hssr.insert(0, "1.4")
        self.hssr.grid(row=3, column=2, sticky=W)

        Label(master, text="  ").grid(row=4, columnspan=5, sticky=N)

        self.useScale = Tkinter.IntVar()
        Tkinter.Checkbutton(master, text="Scale HS %sG by the Minimum Distance to the Ligand" % delta, variable=self.useScale, onvalue=1, offvalue=0).grid(row=5, columnspan=5, sticky=W)

        Label(master, text="  ").grid(row=10, columnspan=5, sticky=N)

    def read_base_dir(self):
        base_dir = ""
        base_dir = askdirectory(title="Select base directory", mustexist=1)
        self.eproj.delete(0, END)
        self.eproj.insert(0, base_dir)
    def apply(self):
        self.project_dir = "%s" % (self.eproj.get())
        self.HSSradius = float(self.hssr.get())
        self.useScale = self.useScale.get()
        self.ok_flag = 1

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

class displayResultsDialog:
    def __init__(self, parent):
        top = self.top = Toplevel(parent)
        self.ok_flag = 0
                
        w1 = Tkinter.Label(top,
                          text = MAL_INFO, 
                          background = 'black',
                          foreground = 'white')
        w1.grid(row=0, columnspan=5, sticky='we')
        w1.grid_columnconfigure(0, weight=1)

        Label(top, text="Protein Desolvation Energy Estimated", foreground='blue', font=('Times', '12', 'bold')).grid(row=1, columnspan=4, sticky=N)
        Label(top, text="(Unit: kcal/mol)", font=('Times', '12')).grid(row=2, columnspan=4, sticky=N)
        Label(top, text="  ").grid(row=3, columnspan=4, sticky=N)
                
        Label(top, text="Ligand name        ", font=("Times", "12")).grid(row=10, column=0, sticky=W)
        Label(top, text="%sG "%delta, font=("Times", "12")).grid(row=10, column=1, sticky=E)
        Label(top, text="-T%sS "%delta, font=("Times", "12")).grid(row=10, column=2, sticky=E)
        Label(top, text="%sH "%delta, font=("Times", "12")).grid(row=10, column=3, sticky=E)

        ci = 0
        for na in range(num_ligands):
            Label(top, text=liglist[na], font=("Times", "12")).grid(row=ci+11, column=0)#, sticky=WE)
            Label(top, text="%.2f " % desol_G[na], font=("Times", "12")).grid(row=ci+11, column=1, sticky=E)
            Label(top, text="%.2f " % desol_S[na], font=("Times", "12")).grid(row=ci+11, column=2, sticky=E)
            Label(top, text="%.2f " % desol_H[na], font=("Times", "12")).grid(row=ci+11, column=3, sticky=E)
            ci += 1
        b = Button(top, text="Close", padx=5, command=self.ok)        
        b.grid(row=num_ligands+30, column=1, sticky=S)

    def ok(self):
        self.top.destroy()



read_settings()

