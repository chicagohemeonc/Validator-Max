#!/usr/bin/env python

import unittest
from PeptideFragmentSingleton import PeptideFragment

f = PeptideFragment( 1 ) # mass type = mono

class TestPeptideFragmentationWithMods(unittest.TestCase):

    def setUp(self):
        #                  charge,  modifications
        f.analyze('ACKRM', 3,       "140@1 16.0@M 57.1@3")

    def test_sequence(self):
        self.assertEqual(f.sequence(), "ACKRM")

    def test_composition(self):
        self.assertEqual(f.composition(), "C(23) H(43) N(9) O(5) S(2)")

    def test_peptide_mass(self):
        l = [item for item in f.peptide_mass()]
        self.assertEqual(l[0], 820.3934217640001)
        self.assertEqual(l[1], 821.40069823088)
        self.assertEqual(l[2], 411.20398734888)
        self.assertEqual(l[3], 274.47175038821337)
        self.assertEqual(l[4], 206.10563190788)
        self.assertEqual(l[5], 165.08596081968)
        self.assertEqual(l[6], 137.73951342754665)
        self.assertEqual(l[7], 118.20633671888001)
        self.assertEqual(l[8], 103.55645418738)
        self.assertEqual(l[9], 92.16210110732445)
        self.assertEqual(l[10], 83.04661864328)

    def test_a_ions(self):
        l = [item for item in f.a_ions()]
        self.assertEqual(l[0], 62.02134285887999)
        self.assertEqual(l[1], 96.35773769421333)
        self.assertEqual(l[2], 158.08939204501334)
        self.assertEqual(l[3], 210.12309572914663)
        self.assertEqual(l[4], 0.0)

    def test_b_ions(self):
        l = [item for item in f.b_ions()]
        self.assertEqual(l[0], 71.35298106888)
        self.assertEqual(l[1], 105.68937590421332)
        self.assertEqual(l[2], 167.42103025501333)
        self.assertEqual(l[3], 219.45473393914665)
        self.assertEqual(l[4], 0.0)

    def test_c_ions(self):
        l = [item for item in f.c_ions()]
        self.assertEqual(l[0], 77.02849743741332)
        self.assertEqual(l[1], 111.36489227274666)
        self.assertEqual(l[2], 173.0965466235467)
        self.assertEqual(l[3], 225.13025030767997)
        self.assertEqual(l[4], 0.0)

    def test_x_ions(self):
        l = [item for item in f.x_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 212.78580063941334)
        self.assertEqual(l[2], 178.44940580408002)
        self.assertEqual(l[3], 116.71775145328002)
        self.assertEqual(l[4], 64.68404776914669)

    def test_y_ions(self):
        l = [item for item in f.y_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 204.12604578621335)
        self.assertEqual(l[2], 169.78965095088003)
        self.assertEqual(l[3], 108.05799660008002)
        self.assertEqual(l[4], 56.02429291594669)

    def test_z_ions(self):
        l = [item for item in f.z_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 198.45052941768003)
        self.assertEqual(l[2], 164.11413458234668)
        self.assertEqual(l[3], 102.3824802315467)
        self.assertEqual(l[4], 50.34877654741336)

    def test_zdot_ions(self):
        l = [item for item in f.zdot_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 198.78647109608002)
        self.assertEqual(l[2], 164.4500762607467)
        self.assertEqual(l[3], 102.7184219099467)
        self.assertEqual(l[4], 50.684718225813356)


class TestPeptideFragmentationWithNTerminalMod(unittest.TestCase):

    def setUp(self):
        #                  charge,  modifications
        f.analyze('ACKRM', 3,       "16.0@[")

    def test_sequence(self):
        self.assertEqual(f.sequence(), "ACKRM")

    def test_composition(self):
        self.assertEqual(f.composition(), "C(23) H(43) N(9) O(5) S(2)")

    def test_peptide_mass(self):
        l = [item for item in f.peptide_mass()]
        self.assertEqual(l[0], 623.2934217640001)
        self.assertEqual(l[1], 624.30069823088)
        self.assertEqual(l[2], 312.65398734888)
        self.assertEqual(l[3], 208.77175038821335)
        self.assertEqual(l[4], 156.83063190788)
        self.assertEqual(l[5], 125.66596081968001)
        self.assertEqual(l[6], 104.88951342754666)
        self.assertEqual(l[7], 90.04919386173715)
        self.assertEqual(l[8], 78.91895418738)
        self.assertEqual(l[9], 70.26210110732444)
        self.assertEqual(l[10], 63.33661864328)

    def test_a_ions(self):
        l = [item for item in f.a_ions()]
        self.assertEqual(l[0], 20.688009525546665)
        self.assertEqual(l[1], 55.02440436087999)
        self.assertEqual(l[2], 97.72272537834668)
        self.assertEqual(l[3], 149.75642906248)
        self.assertEqual(l[4], 0.0)

    def test_b_ions(self):
        l = [item for item in f.b_ions()]
        self.assertEqual(l[0], 30.019647735546666)
        self.assertEqual(l[1], 64.35604257087999)
        self.assertEqual(l[2], 107.05436358834667)
        self.assertEqual(l[3], 159.08806727248)
        self.assertEqual(l[4], 0.0)

    def test_c_ions(self):
        l = [item for item in f.c_ions()]
        self.assertEqual(l[0], 35.69516410408)
        self.assertEqual(l[1], 70.03155893941333)
        self.assertEqual(l[2], 112.72987995688)
        self.assertEqual(l[3], 164.76358364101335)
        self.assertEqual(l[4], 0.0)

    def test_x_ions(self):
        l = [item for item in f.x_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 188.41913397274666)
        self.assertEqual(l[2], 154.08273913741334)
        self.assertEqual(l[3], 111.38441811994669)
        self.assertEqual(l[4], 59.35071443581336)

    def test_y_ions(self):
        l = [item for item in f.y_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 179.75937911954668)
        self.assertEqual(l[2], 145.42298428421336)
        self.assertEqual(l[3], 102.72466326674669)
        self.assertEqual(l[4], 50.69095958261335)

    def test_z_ions(self):
        l = [item for item in f.z_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 174.08386275101336)
        self.assertEqual(l[2], 139.74746791568)
        self.assertEqual(l[3], 97.04914689821335)
        self.assertEqual(l[4], 45.01544321408002)

    def test_zdot_ions(self):
        l = [item for item in f.zdot_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 174.41980442941335)
        self.assertEqual(l[2], 140.08340959408002)
        self.assertEqual(l[3], 97.38508857661337)
        self.assertEqual(l[4], 45.35138489248002)


class TestPeptideFragmentationWithCTerminalMod(unittest.TestCase):

    def setUp(self):
        #                  charge,  modifications
        f.analyze('ACKRM', 3,       "20.0@]")

    def test_sequence(self):
        self.assertEqual(f.sequence(), "ACKRM")

    def test_composition(self):
        self.assertEqual(f.composition(), "C(23) H(43) N(9) O(5) S(2)")

    def test_peptide_mass(self):
        l = [item for item in f.peptide_mass()]
        self.assertEqual(l[0], 627.2934217640001)
        self.assertEqual(l[1], 628.30069823088)
        self.assertEqual(l[2], 314.65398734888)
        self.assertEqual(l[3], 210.10508372154666)
        self.assertEqual(l[4], 157.83063190788)
        self.assertEqual(l[5], 126.46596081968)
        self.assertEqual(l[6], 105.55618009421333)
        self.assertEqual(l[7], 90.62062243316572)
        self.assertEqual(l[8], 79.41895418738)
        self.assertEqual(l[9], 70.7065455517689)
        self.assertEqual(l[10], 63.736618643279996)

    def test_a_ions(self):
        l = [item for item in f.a_ions()]
        self.assertEqual(l[0], 15.354676192213333)
        self.assertEqual(l[1], 49.691071027546656)
        self.assertEqual(l[2], 92.38939204501334)
        self.assertEqual(l[3], 144.4230957291467)
        self.assertEqual(l[4], 0.0)

    def test_b_ions(self):
        l = [item for item in f.b_ions()]
        self.assertEqual(l[0], 24.68631440221333)
        self.assertEqual(l[1], 59.02270923754666)
        self.assertEqual(l[2], 101.72103025501333)
        self.assertEqual(l[3], 153.7547339391467)
        self.assertEqual(l[4], 0.0)

    def test_c_ions(self):
        l = [item for item in f.c_ions()]
        self.assertEqual(l[0], 30.361830770746664)
        self.assertEqual(l[1], 64.69822560607999)
        self.assertEqual(l[2], 107.39654662354667)
        self.assertEqual(l[3], 159.43025030768)
        self.assertEqual(l[4], 0.0)

    def test_x_ions(self):
        l = [item for item in f.x_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 195.08580063941335)
        self.assertEqual(l[2], 160.74940580408)
        self.assertEqual(l[3], 118.05108478661334)
        self.assertEqual(l[4], 66.01738110248003)

    def test_y_ions(self):
        l = [item for item in f.y_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 186.42604578621334)
        self.assertEqual(l[2], 152.08965095088)
        self.assertEqual(l[3], 109.39132993341336)
        self.assertEqual(l[4], 57.357626249280024)

    def test_z_ions(self):
        l = [item for item in f.z_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 180.75052941768)
        self.assertEqual(l[2], 146.4141345823467)
        self.assertEqual(l[3], 103.71581356488002)
        self.assertEqual(l[4], 51.68210988074669)

    def test_zdot_ions(self):
        l = [item for item in f.zdot_ions()]
        self.assertEqual(l[0], 0.0)
        self.assertEqual(l[1], 181.08647109608)
        self.assertEqual(l[2], 146.75007626074668)
        self.assertEqual(l[3], 104.05175524328003)
        self.assertEqual(l[4], 52.01805155914669)

if __name__ == '__main__':
    unittest.main()
