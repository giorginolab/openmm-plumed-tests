from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

from openmmplumed import PlumedForce


prmtop = AmberPrmtopFile('w.prmtop')
inpcrd = AmberInpcrdFile('w.inpcrd')
#pdb = PDBFile("w.pdb")




system = prmtop.createSystem(nonbondedMethod=PME, 
	nonbondedCutoff=9*angstrom,  constraints=HBonds)



#r: RESTRAINT ARG=d AT=5 KAPPA=150
#UNITS LENGTH=A 

plumedScript = """
d: DISTANCE ATOMS=1,2758
r: RESTRAINT ARG=d AT=.5 KAPPA=1500
PRINT ARG=* FILE=COLVAR
"""

system.addForce(PlumedForce(plumedScript))

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

req_plt = Platform.getPlatformByName('CUDA')

simulation = Simulation(prmtop.topology, system, integrator, req_plt)
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

simulation.saveState("pre.xml")

print("Minimizing")
simulation.minimizeEnergy(maxIterations=100)
simulation.saveState("minimized.xml")

print("Running")
simulation.reporters.append(DCDReporter('output.dcd', 100))
simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
                    potentialEnergy=True, temperature=True))
simulation.step(10000)

simulation.saveState("post.xml")

