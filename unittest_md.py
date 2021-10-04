#!/usr/bin/env python3

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
import sys, unittest
from md import calcenergy
# Use Asap for a huge performance increase if it is installed
use_asap = True

if use_asap:
    from asap3 import EMT
    size = 10
else:
    from ase.calculators.emt import EMT
    size = 3

class MdTests(unittest.TestCase):

    def test_calcenergy(self):
        interval=10
        size = 10
        # Set up a crystal
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(size, size, size),
                                  pbc=True)

        MaxwellBoltzmannDistribution(atoms, temperature_K=300)
        atoms.calc = EMT()
        energies = calcenergy(atoms)
        self.assertAlmostEqual(300, energies[2], delta=10)
        dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

        dyn.run(200)
        
        energies = calcenergy(atoms)
        #check the potential energy
        self.assertAlmostEqual(0.020, energies[0], delta=0.01)
        # check the kinetic energy
        self.assertAlmostEqual(0.020, energies[1], delta=0.01)
        #self.assertEqual(300, energies[2])
        #self.assertTrue(True)

if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(MdTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
