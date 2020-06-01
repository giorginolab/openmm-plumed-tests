
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

from openmmplumed import PlumedForce

import numpy as np


N = 4

system=System()

positions = []
for i in range(N):
    system.addParticle(1.0)
    positions.append( [i, 0.1*i, -.3*i] )
positions=np.array(positions)

script = """
d: DISTANCE ATOMS=1,3
BIASVALUE ARG=d
"""

plumed = PlumedForce(script)
system.addForce(plumed)

integ=LangevinIntegrator(300., 1., 1.,)

platform = Platform.getPlatformByName("CUDA")
context = Context(system, integ, platform)
context.setPositions(positions)


state = context.getState(getForces=True, getEnergy=True)

delta = positions[0, :]-positions[2, :]

dist = np.sqrt(delta.dot(delta))

print(dist, state.getPotentialEnergy())
print(-delta/dist, state.getForces())
