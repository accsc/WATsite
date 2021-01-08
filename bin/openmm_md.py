#!/usr/bin/env python
#title           :openmm_md.py
#description     :This will perform MD simulation with OpenMM package
#author          :Ying Yang
#date            :11-13-2018
#usage           :python openmm_md.py -h
#python_version  :2.7  
#==============================================================================

import parmed as pmd
from parmed import unit as u
import mdtraj
from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app.internal.unitcell import computeLengthsAndAngles
from collections import defaultdict, OrderedDict

from multiprocessing import Process
import threading

from argparse import ArgumentParser
import numpy as np
import os, sys, time, logging, math
t0 = time.time()
logging.basicConfig()


#================================#
#       Set up the logger        #
#================================#

logger = logging.getLogger(__name__)
logger.setLevel('INFO')
sh = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(fmt='%(asctime)s - %(message)s', datefmt="%H:%M:%S")
sh.setFormatter(formatter)
logger.addHandler(sh)
logger.propagate = False

#================================#
#          Subroutines           #
#================================#

'''class EnergyReporter(object):
  def __init__(self, fileout, prev_frame, reportInterval, waterlist, watermodel):
    self._reportInterval = reportInterval
    self._waterlist = waterlist
    self._watermodel = watermodel
    self._fileout = fileout
    self._prev_frame = prev_frame
  #def __del__(self):
  #  self._out.close()
  #  return
  def describeNextReport(self, simulation):
    steps = self._reportInterval - simulation.currentStep%self._reportInterval
    return (steps, False, False, False, True) 
  def report(self, simulation, state):
    epp = simulation.context.getState(getEnergy=True).getPotentialEnergy_PerParticle()
    #oneprocess = Process(target=writeEnergyPerParticle, args=(epp.value_in_unit(u.kilocalories_per_mole), self._fileout, self._prev_frame, simulation.currentStep, self._reportInterval, self._watermodel, self._waterlist))
    #threads.append(oneprocess)    
    #oneprocess.start() # will run 
    return epp

def writeEnergyPerParticle(epp,filename,prev_frame,simulation_step,interval,watermodel,waterlist):
  frameNum = int(prev_frame) + int(simulation_step) / int(interval)
  write2file= open(filename +'_'+ str(frameNum), "w")
  write2file.write("$STEP: % 12d\n" % simulation_step)
  site3_water = ["SPC/E", "TIP3P", "AMOEBA"]
  site4_water = ["OPC", "TIP4P", "TIP4PEW"]
  site5_water = ["TIP5P"]
  if watermodel in site3_water:
    numSite = 3
  elif watermodel in site4_water:
    numSite = 4
  elif watermodel in site5_water:
    numSite = 5
  else:
    print "Error: Water model '%s' not supported..." % watermodel
  for i in range(len(epp)):
    j = i+1
    if j in waterlist:
      sumen = 0.0;
      for k in range(numSite):
        sumen += epp[i+k]
      write2file.write("% 10d % 12.5f\n" % (j, sumen))  
  write2file.close()    
  return 0'''

def createListOfWaterOxygens(pdbfile):
  waterlist = []
  filein = open(pdbfile)
  for line in filein:
    if line.find("ATOM") == 0 or line.find("HETATM") == 0:
      if line[12:20].find("O   WAT") >= 0 or line[12:20].find("O   HOH") >= 0  or line[12:20].find("OW   HOH")>=0 or line[12:20].find("OW  SOL") >= 0 :
        atomid = int(line[6:11])
        waterlist.append(atomid)
  return waterlist

def MTSVVVRIntegrator(temperature, collision_rate, timestep, system, ninnersteps=4):
    """
    Create a multiple timestep velocity verlet with velocity randomization (VVVR) integrator.
    
    ARGUMENTS

    temperature (numpy.unit.Quantity compatible with kelvin) - the temperature
    collision_rate (numpy.unit.Quantity compatible with 1/picoseconds) - the collision rate
    timestep (numpy.unit.Quantity compatible with femtoseconds) - the integration timestep
    system (simtk.openmm.System) - system whose forces will be partitioned
    ninnersteps (int) - number of inner timesteps (default: 4)

    RETURNS

    integrator (openmm.CustomIntegrator) - a VVVR integrator

    NOTES
    
    This integrator is equivalent to a Langevin integrator in the velocity Verlet discretization with a
    timestep correction to ensure that the field-free diffusion constant is timestep invariant.  The inner
    velocity Verlet discretization is transformed into a multiple timestep algorithm.

    REFERENCES

    VVVR Langevin integrator: 
    * http://arxiv.org/abs/1301.3800
    * http://arxiv.org/abs/1107.2967 (to appear in PRX 2013)    
    
    TODO

    Move initialization of 'sigma' to setting the per-particle variables.
    
    """
    # Multiple timestep Langevin integrator.
    for i in system.getForces():
        if i.__class__.__name__ in ["NonbondedForce", "CustomNonbondedForce", "AmoebaVdwForce", "AmoebaMultipoleForce"]:
            # Slow force.
            print i.__class__.__name__, "is a Slow Force"
            i.setForceGroup(1)
        else:
            print i.__class__.__name__, "is a Fast Force"
            # Fast force.
            i.setForceGroup(0)

    kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
    kT = kB * temperature
    
    integrator = mm.CustomIntegrator(timestep)
    
    integrator.addGlobalVariable("dt_fast", timestep/float(ninnersteps)) # fast inner timestep
    integrator.addGlobalVariable("kT", kT) # thermal energy
    integrator.addGlobalVariable("a", np.exp(-collision_rate*timestep)) # velocity mixing parameter
    integrator.addGlobalVariable("b", np.sqrt((2/(collision_rate*timestep)) * np.tanh(collision_rate*timestep/2))) # timestep correction parameter
    integrator.addPerDofVariable("sigma", 0) 
    integrator.addPerDofVariable("x1", 0) # position before application of constraints

    #
    # Pre-computation.
    # This only needs to be done once, but it needs to be done for each degree of freedom.
    # Could move this to initialization?
    #
    integrator.addComputePerDof("sigma", "sqrt(kT/m)")

    # 
    # Velocity perturbation.
    #
    integrator.addComputePerDof("v", "sqrt(a)*v + sqrt(1-a)*sigma*gaussian")
    integrator.addConstrainVelocities();
    
    #
    # Symplectic inner multiple timestep.
    #
    integrator.addUpdateContextState(); 
    integrator.addComputePerDof("v", "v + 0.5*b*dt*f1/m")
    for innerstep in range(ninnersteps):
        # Fast inner symplectic timestep.
        integrator.addComputePerDof("v", "v + 0.5*b*dt_fast*f0/m")
        integrator.addComputePerDof("x", "x + v*b*dt_fast")
        integrator.addComputePerDof("x1", "x")
        integrator.addConstrainPositions();        
        integrator.addComputePerDof("v", "v + 0.5*b*dt_fast*f0/m + (x-x1)/dt_fast")
    integrator.addComputePerDof("v", "v + 0.5*b*dt*f1/m") # TODO: Additional velocity constraint correction?
    integrator.addConstrainVelocities();

    #
    # Velocity randomization
    #
    integrator.addComputePerDof("v", "sqrt(a)*v + sqrt(1-a)*sigma*gaussian")
    integrator.addConstrainVelocities();

    return integrator


parser = ArgumentParser()
#\x1b[1;91m(Required)\x1b[0m

group = parser.add_argument_group('Input File Options')
group.add_argument('-p', '--pdb', dest='pdb',  metavar='<PDB FILE>', help='PDB file with target system.')
group.add_argument('-i', '--inpcrd', dest='crd', metavar='<CRD FILE>', help='Amber inpcrd file with target system.')
group.add_argument('-t', '--topology', dest='top', metavar='<TOP FILE>', help='Amber prmtop file with target system')
group.add_argument('--xml', dest='xml', default='system.xml', metavar='FILE', help='''OpenMM System XML file. Default is %(default)s''')
group.add_argument('-s', '--state', dest='state', metavar='FILE', required=True, help='''Restart file (format: chk, xml, rst)\x1b[1;91m(Required)\x1b[0m''')
group.add_argument('-c', '--cuda', dest='cuda', default='0', metavar='INT', help='''Index of CUDA device: 0, 1, 2...''')
group.add_argument('-w', '--water', dest='wat_model', help='Water model: AMOEBA; SPC/E; TIP3P; OPC; TIP4P; TIP4PEW.')

group = parser.add_argument_group('Positional Restraint / Repulsion for Probe Options')
group.add_argument('--reference', dest='reference', metavar='PDB FILE', help='restraint reference PDB structure (default None)', default=None)
group.add_argument('--restrain-mask', dest='restraints', metavar='MASK', help='restraint mask (default None)', default=None)
group.add_argument('-k', '--force-constant', dest='force_constant', type=float, metavar='FLOAT', help='''Force constant for cartesian constraints. Default 10 kcal/mol/A^2''', default=10)
group.add_argument('--repulsion-mask', dest='repulsion', metavar='MASK', help='amber atom mask for probe repulsion. E.g., :BBZ@C1,C2,C3,C4,C5,C6s (default None)', default=None)

group = parser.add_argument_group('Output File Options')
group.add_argument('-o' , '--output', dest='output', default=sys.stdout, metavar='FILE', help='''Output file for energies''')
group.add_argument('-r' , '--report', dest='rep', default='sys_NVT.info', metavar='FILE', help='''Output file for progress''')
group.add_argument('-x', '--trajectory', dest='trajectory', default='sys_NVT.nc', metavar='FILE', help='''NetCDF trajectory to generate. Snapshots written every --interval steps.''')
group.add_argument('--restart', dest='restart', default='sys_NVT.xml', metavar='FILE', help='''Xml serialized format with information to restart the simulation with another run''')
group.add_argument('-e', '--wat_egy', dest='watEgyPrefix', default='waterEnergy', metavar='FILE', help='Output file for water energies. Default is %(default)s')
group.add_argument('--prev_frame', dest='prev_frame', default=0, metavar='INT', type=int, help='Last frame number from previous run. Default is %(default)s')


group = parser.add_argument_group('Simulation Options')
group.add_argument('-n', '--num-steps', dest='num_steps', required=True, type=int, help='Number of MD steps to run. \x1b[35m(Required)\x1b[0m', metavar='INT')
group.add_argument('--interval', dest='interval', default=500, metavar='INT', help='Interval between printing state data. Default 500', type=int)
group.add_argument('--dt', dest='timestep', type=float, metavar='FLOAT', help='''time step for integrator (outer time-step for RESPA integrator) Default 1 fs''', default=1.0)
group.add_argument('--temp', dest='temp', type=float, metavar='FLOAT', help='''target temperature for NVT simulation. Default %(default)s K''', default=300.0)
group.add_argument('--temp1', dest='temp1', type=int, metavar='FLOAT', help='''initial temperature for annealing.''')
group.add_argument('--temp2', dest='temp2', type=int, metavar='FLOAT', help='''final   temperature for annealing.''')
group.add_argument('--npt', dest='npt', default=False, action='store_true', help='Do constant pressure simulation')
group.add_argument('--random', dest='random', default=False, action='store_true', help='Initial velocity by random seed')


group = parser.add_argument_group('Integrator Options')
group.add_argument('--gamma_ln', dest='gamma_ln', type=float, default=1.0, metavar='FLOAT', help='''collision frequency for \x1b[1;91mLangevinIntegrator\x1b[0m. Default %(default)s ps-1''')
group.add_argument('--nrespa', dest='nrespa', type=int, metavar='INT', help='''Number of inner time steps to run (evaluating fast forces) for every outer timestep. For \x1b[1;91mMTS or MTSVVVR\x1b[0m integrator.''')
group.add_argument('--mtsvvvr', dest='mtsvvvr', default=False, action='store_true', help='Use \x1b[1;91mMTSVVVRIntegrator\x1b[0m adapted from Lee-Ping Wang')

########################################################################

opt = parser.parse_args()

logger.info('Command line:\n\t%s' % ' '.join(sys.argv))
logger.info('Parsing system XML file %s' % opt.xml)
with open(opt.xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
    
#for i in range(system.getNumForces()):
#    logger.info('Force # %i -- %s' % (i, system.getForce(i)) )

if opt.crd and opt.top:
    prmtop = app.AmberPrmtopFile(opt.top)
    inpcrd = app.AmberInpcrdFile(opt.crd)
    pos = inpcrd.positions
    top = prmtop.topology

elif opt.pdb:
    pdb = app.PDBFile(opt.pdb)
    #pdb = pmd.load_file(opt.pdb)
    pos = pdb.positions
    top = pdb.topology


########################################################################
if opt.restraints or opt.repulsion:
    ref_pdb = pmd.load_file(opt.reference)

# Add cartesian restraints if desired
if opt.restraints:

    print('Adding restraints (k=%s kcal/mol/A^2) from %s' % (opt.force_constant, opt.restraints))
    sel = pmd.amber.AmberMask(ref_pdb, opt.restraints).Selection()
    const = opt.force_constant * u.kilocalories_per_mole/u.angstroms**2
    const = const.value_in_unit_system(u.md_unit_system)
    force = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    force.addGlobalParameter('k', const)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    for i, atom_crd in enumerate(ref_pdb.positions):
        if sel[i]:
            force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)

# Add repulsion force between probes if desired
if opt.repulsion:
    excludedIxn = defaultdict(list)
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if force.__class__.__name__ == 'NonbondedForce':
            for j in range(force.getNumExceptions()):
                (p1, p2, q, sig, eps) = force.getExceptionParameters(j)
                excludedIxn[p1].append(p2)
    #print(excludedIxn.items())

    const = 2.5 * u.kilocalories_per_mole/u.angstroms**2  # Need to test K and min_distance
    const = const.value_in_unit_system(u.md_unit_system)
    min_dis = 6.8 * u.angstroms
    min_dis = min_dis.value_in_unit_system(u.md_unit_system)

    flat_bottom_force = mm.CustomBondForce('step(r0-r) * k * (r-r0)^2')
    flat_bottom_force.addGlobalParameter('k', const)
    flat_bottom_force.addGlobalParameter('r0', min_dis)

    repulsiveAtoms = []
    sel_atom = pmd.amber.AmberMask(ref_pdb, opt.repulsion).Selection() # ":F29@C1,C2,C3,C4,C5,C6" return a list 0 or 1. Selected atoms are 1. 
    for i, atom in enumerate(ref_pdb.atoms):
        if sel_atom[i]:
            #print("i=",i, "    ", atom)
            repulsiveAtoms.append(i)
    for i in repulsiveAtoms:
        for j in repulsiveAtoms:
            if i != j and not (j in excludedIxn[i]):
                flat_bottom_force.addBond(i, j)
    system.addForce(flat_bottom_force)

########################################################################
if opt.npt:
### Perform Pressure coupling
    logger.info('Using isotropic barostat')
    baro = mm.MonteCarloBarostat(1*u.bar, opt.temp*u.kelvin)
    system.addForce(baro)

### Set integrator 
if opt.nrespa:
    if opt.mtsvvvr:
        logger.info("Use MTSVVVRIntegrator: a multiple timestep Langevin integrator with %.2f / %.2f fs outer/inner timestep." % (opt.timestep, opt.timestep/opt.nrespa ) )
        logger.info("The stochastic thermostat collision frequency is %s ps^-1" % opt.gamma_ln )
        integrator = MTSVVVRIntegrator(opt.temp*u.kelvin, opt.gamma_ln/u.picosecond, opt.timestep*u.femtosecond, system, opt.nrespa)

    else:
        logger.info("Use RESPA MTS Integrator with Anderson Thermostat. %.2f / %.2f fs outer/inner timestep." % (opt.timestep, opt.timestep/opt.nrespa ) )
        logger.info('Adding Anderson thermostat at %s K, collison rate of %s psec^-1' % (opt.temp, opt.gamma_ln) ) 
        thermo = mm.AndersenThermostat(opt.temp*u.kelvin, opt.gamma_ln/u.picosecond)
        system.addForce(thermo)

        slow = (mm.AmoebaMultipoleForce, mm.AmoebaVdwForce, mm.AmoebaGeneralizedKirkwoodForce, mm.AmoebaWcaDispersionForce, mm.NonbondedForce, mm.CustomNonbondedForce )
        found_slow = False
        for force in system.getForces():
            if isinstance(force, slow):
                found_slow = True
                force.setForceGroup(1)
            else:
                force.setForceGroup(0)
        if not found_slow:
            raise ValueError('No slow forces found for MTS integrator!')
        integrator = mm.MTSIntegrator(opt.timestep*u.femtoseconds, [(0, opt.nrespa), (1, 1)])
else:
    ### Temperature coupling with Langevin integrator
    integrator = mm.LangevinIntegrator(opt.temp*u.kelvin, opt.gamma_ln/u.picosecond, opt.timestep*u.femtoseconds)
    logger.info('Use Langevin integrator: %8.2fK, %8.2f ps-1, %8.2f fs' % (opt.temp, opt.gamma_ln, opt.timestep) )

#for i in range(system.getNumForces()):
#    logger.info('Force # %i -- %s' % (i, system.getForce(i)) )

########################################################################
platform = mm.Platform_getPlatformByName("CUDA")
platformProperties = {}
platformProperties['CudaPrecision'] = 'single'
platformProperties['CudaDeviceIndex'] = opt.cuda

sim = app.Simulation(top, system, integrator, platform, platformProperties )

sim.reporters.append(
        pmd.openmm.StateDataReporter(opt.output, reportInterval=opt.interval, volume=True,density=True,separator='\t') )
sim.reporters.append(
        pmd.openmm.ProgressReporter(opt.rep, opt.interval, opt.num_steps) )
sim.reporters.append(
        mdtraj.reporters.NetCDFReporter(opt.trajectory, opt.interval) )

#sim.reporters.append(app.DCDReporter(opt.trajectory, opt.interval))


if opt.wat_model:

    #system.set_
    #print system.waterlist
    site3_water = ["SPC/E", "TIP3P", "AMOEBA"]
    site4_water = ["OPC", "TIP4P", "TIP4PEW"]
    site5_water = ["TIP5P"]
    if opt.wat_model in site3_water:
      numSite = 3
    elif opt.wat_model in site4_water:
      numSite = 4
    elif opt.wat_model in site5_water:
      numSite = 5
      
    system.set_watermodel(numSite)
    waterlist = createListOfWaterOxygens(opt.pdb)
    for i in waterlist:
      system.add_waterlist(i)
    system.set_intervalwater(opt.interval)
    system.set_prevframe(opt.prev_frame)

    logger.info("Read water index from file: %s" % (opt.pdb))
    logger.info("Water energy will be written to files with prefix: %s" % opt.watEgyPrefix)

    watEgyFolder = "waterEnergy"
    try:
        os.mkdir(watEgyFolder)
    except OSError as e:
        if e.errno == 17: # errno.EEXIST
            logger.info("Warning: (%s) exists" % watEgyFolder)
        else:
            raise
    #threads = []
    
else:
    system.set_watermodel(0)
    system.set_intervalwater(opt.interval)


#amr is messing in watsite

########################################################################
# Print out some more information about the system
logger.info("--== System Information ==--")
logger.info("Number of particles   : %i" % sim.context.getSystem().getNumParticles())
logger.info("Number of constraints : %i" % sim.context.getSystem().getNumConstraints())
for f in sim.context.getSystem().getForces():
    if f.__class__.__name__ == 'AmoebaMultipoleForce':
        f.updateParametersInContext(sim.context)
        logger.info("AMOEBA PME order      : %i" % f.getPmeBSplineOrder())
        logger.info("AMOEBA PME grid       : %s" % str(f.getPmeGridDimensions()))
    if f.__class__.__name__ == 'NonbondedForce':
        method_names = ["NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME"]
        logger.info("Nonbonded method      : %s" % method_names[f.getNonbondedMethod()])
        logger.info("Number of particles   : %i" % f.getNumParticles())
        logger.info("Number of exceptions  : %i" % f.getNumExceptions())
        if f.getNonbondedMethod() > 0:
            logger.info("Nonbonded cutoff      : %.3f nm" % (f.getCutoffDistance() / u.nanometer))
            if f.getNonbondedMethod() >= 3:
                logger.info("Ewald error tolerance : %.3e" % (f.getEwaldErrorTolerance()))
            logger.info("LJ switching function : %i" % f.getUseSwitchingFunction())
            if f.getUseSwitchingFunction():
                logger.info("LJ switching distance : %.3f nm" % (f.getSwitchingDistance() / u.nanometer))

########################################################################

# Read the coordinates & verlocities from the restart file
if opt.state is not None:
    logger.info('Setting coordinates and velocities from restart file %s' % opt.state)
    if opt.state[-3:] == 'xml':
        with open(opt.state, 'r') as f:
            sim.context.setState(mm.XmlSerializer.deserialize(f.read()))
    elif opt.state[-3:] == 'chk':
        sim.loadCheckpoint(opt.state)
    else:
    #jason's code that is supposed to work for any restart file type:
        rst = pmd.load_file(opt.state)
        sim.context.setPositions(rst.coordinates[-1]*u.angstroms)
        sim.context.setVelocities(rst.velocities[-1]*u.angstroms/u.picoseconds)
        sim.context.setPeriodicBoxVectors(*pmd.geometry.box_lengths_and_angles_to_vectors(*rst.box))
        if hasattr(rst, 'time'):
            try:
                sim.context.setTime(rst.time[-1])
            except TypeError:
                sim.context.setTime(rst.time)
else:
    logger.info('Failed to restart from minimized conformation...\nCheck input!\n')
    sys.exit(0)
sim.context.applyConstraints(1e-5)


########################################################################

# Randomize initial velocity 
if opt.random:
    sim.context.setVelocitiesToTemperature(300*u.kelvin)

########################################################################

# Start the simulation run
t1 = time.time()

if opt.temp1 and opt.temp2:
    perTstep = int(opt.num_steps / abs(opt.temp2-opt.temp1))

    logger.info('\nRunning the annealing simulation from %d k to %d k for %d steps at each temperature!\nTotal number of steps: %d' % (opt.temp1, opt.temp2, perTstep, opt.num_steps) )

    if opt.temp2 > opt.temp1:
        for i in range(opt.temp1, opt.temp2):
            integrator.setTemperature(i*u.kelvin)
            sim.step( perTstep )
    else:
        for i in range(opt.temp1, opt.temp2, -1):
            integrator.setTemperature(i*u.kelvin)
            sim.step( perTstep )
else:
    logger.info('Running the simulation for %d steps!' % opt.num_steps)
    sim.step(opt.num_steps)


########################################################################

# Save the final coordinates into pdb, and write restart file in XML format
print('Finished. Writing serialized XML restart file...')
final_state = sim.context.getState(getPositions=True, getVelocities=True, getForces=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions(), getEnergy=True)
positions = final_state.getPositions()
pbc_box = final_state.getPeriodicBoxVectors()
a, b, c, alpha, beta, gamma = computeLengthsAndAngles(pbc_box)
RAD_TO_DEG = 180/math.pi

outxml = opt.restart
outpdb = opt.restart.split(".")[0]+".pdb"

with open(outxml, 'w') as f:
    f.write( mm.XmlSerializer.serialize(final_state) )

app.PDBFile.writeModel(top, positions, open(outpdb, 'w'))
Fo = open(outpdb, 'a')
Fo.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n" % (a*10, b*10, c*10, alpha*RAD_TO_DEG, beta*RAD_TO_DEG, gamma*RAD_TO_DEG) )
Fo.write("END\n\n")


#=============================================#
#| Run analysis and save restart information |#
#=============================================#
prodtime = time.time() - t1

logger.info('Getting statistics for the production run.')
print "Total wall time: % .4f seconds" % (time.time() - t0)
print "Production wall time: % .4f seconds" % (prodtime)
print "Simulation speed: % .6f ns/day" % (86400*opt.num_steps*opt.timestep*u.femtosecond/u.nanosecond/(prodtime))

#=============================================#
#|     Wait until all processes are done.    |#
#=============================================#

#if opt.wat_model:
#    for thread in threads:
#        thread.join()
