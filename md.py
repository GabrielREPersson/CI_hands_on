#!/usr/bin/env python3

"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import Trajectory

# Use Asap for a huge performance increase if it is installed
use_asap = True

if use_asap:
    from asap3 import EMT
    size = 10
else:
    from ase.calculators.emt import EMT
    size = 3

def calcenergy(atoms):
    """Function to print the potential, kinetic and total energy."""
    epot = atoms.get_potential_energy() / len(atoms)
    ekin = atoms.get_kinetic_energy() / len(atoms)
    insta_temp = ekin / (1.5 * units.kB)
    etot = epot + ekin
    return [epot, ekin, insta_temp, etot]


def run_md():
    interval=10

    # Set up a crystal
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol="Cu",
                              size=(size, size, size),
                              pbc=True)

    # Describe the interatomic interactions with the Effective Medium Theory
    atoms.calc = EMT()

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
    traj = Trajectory("cu.traj", "w", atoms)
    dyn.attach(traj.write, interval=interval)

    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        energys = calcenergy(atoms)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (energys[0], energys[1], energys[2], energys[3]))


    # Now run the dynamics
    dyn.attach(printenergy, interval=interval)
    printenergy()
    dyn.run(200)

if __name__ == "__main__":
    run_md()
