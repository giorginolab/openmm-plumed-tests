from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

from openmmplumed import PlumedForce

psf = CharmmPsfFile('dia2.psf')
pdb = PDBFile('dia2.pdb')

params = CharmmParameterSet('par_all27_prot_lipid.prm', permissive=True)
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
                          nonbondedCutoff=1*nanometer, constraints=None)

plumedScript = "diala.plumed.nocont"
with open(plumedScript) as f:
    script = f.read()
plumedForce = PlumedForce(script)

req_plt = Platform.getPlatformByName('CUDA')

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)
simulation = Simulation(psf.topology, system, integrator,
                        req_plt, {'DeviceIndex': '1'} )

ctx = simulation.context
platform = ctx.getPlatform()
print(f"Using platform {platform.getName()} with properties:")
for prop in platform.getPropertyNames():
    print(f"    {prop}\t\t{platform.getPropertyValue(ctx,prop)}")

L = 32.0
ctx.setPeriodicBoxVectors([L,0,0], [0,L,0], [0,0,L])

simulation.context.setPositions(pdb.positions)

simulation.saveState("pre.xml")

U=simulation.context.getState(getEnergy=True).getPotentialEnergy()
print("U pre = "+str(U))

simulation.minimizeEnergy(tolerance= 0.001)
minState = simulation.context.getState(getPositions=True, getEnergy=True)

simulation.saveState("minimized.xml")
PDBFile.writeFile(psf.topology, minState.getPositions(),
                  open("minimized.pdb","w"))
U=minState.getPotentialEnergy()
print("U minimized = "+str(U))


# Rebuild the simulation object, this time with Plumed forces.
system.addForce(plumedForce)
integrator2 = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)
simulation2 = Simulation(psf.topology, system, integrator2,
                        req_plt, {'DeviceIndex': '1'} )
simulation2.context.setState(minState)


simulation2.reporters.append(DCDReporter('output.dcd', 1000))
simulation2.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))

simulation2.step(10000000)

