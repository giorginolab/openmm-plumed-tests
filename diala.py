from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

psf = CharmmPsfFile('dia2.psf')
pdb = PDBFile('dia2.pdb')

params = CharmmParameterSet('par_all27_prot_lipid.prm', permissive=True)
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)
simulation = Simulation(psf.topology, system, integrator)

ctx = simulation.context
platform = ctx.getPlatform()
print(f"Using platform {platform.getName()} with properties:")
for prop in platform.getPropertyNames():
    print(f"    {prop}\t\t{platform.getPropertyValue(ctx,prop)}")


simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000000)
