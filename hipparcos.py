import sys, math, pdb

hip_data_path = '/Volumes/Pallas/Hipparcos Star Data/hip_main.dat'

def stars(mag=20, maxdec=math.pi, mindec=-math.pi):
    """A generator-based reader for stars in the Hipparcos catalogue.

    This will return a generator for stars in the Hipparcos catalog.
    Yielded values will be a tuple of the form (hip_id, mag, ra, dec,
    bv), where the hip_id (Hipparcos ID) is an integer, and the
    remaining quantities (Johnson magnitude, RA and Dec in radians,
    and B-V magnitude redness) are floats.
    
    The maxdec and mindec indicate the max and min declination of
    stars that should be returned, respectively.  The declination
    quantities are in radians.  (Declination of pi indicates the
    celestial north pole.)  The default arguments will ensure all
    stars are returned.  Stars are returned in the order read from the
    catalog.  Note that the catalog itself does not contain the stars
    ordered by ID number."""
    f = file(hip_data_path)
    # The fields of the Hipparcos hip_main.dat file are described in page
    # 136 of the document "Contents of the Hipparcos Catalogue," which
    # corresponds to page 180 in the PDF.
    for line in f:
        hip_number = int(line[2:14])
        # There is one troublesome case we must ignore...
        if hip_number == 120412: continue
        johnson_magnitude = float(line[41:46])
        if johnson_magnitude>mag: continue
        ra = int(line[17:19])+int(line[20:22])/60.+float(line[23:28])/3600.
        ra = math.radians(ra * 15)
        dec = int(line[30:32])+int(line[33:35])/60.+float(line[36:40])/3600.
        dec = math.radians(dec)
        if line[29]=='-': dec = -dec
        if dec>maxdec or dec<mindec: continue
        try: bv = float(line[245:251])
        except ValueError: bv = None
        yield hip_number, johnson_magnitude, ra, dec, bv
    f.close()

if __name__ == "__main__":
    import unittest
    class TestStar(unittest.TestCase):
        def testFirst(self):
            """Tests whether the first few stars are correct."""
            s = [(1, 9.10, 1.5999e-05, 0.019007, 0.482),
                 (2, 9.27, 6.6177e-05, -0.34032, 0.999),
                 (3, 6.61, 8.7266e-05,  0.67822,-0.019),
                 (4, 8.06, 1.4617e-04, -0.90571, 0.370),
                 (5, 8.55, 1.7381e-04, -0.70845, 0.902)]
            sgen = stars()
            for a in s:
                b = sgen.next()
                self.assertEqual(len(b), len(a))
                c = zip(a, b)
                for i,j in c:
                    self.assertAlmostEqual(i, j, 5)

        def testCount(self):
            """Test that all HIP IDs are unique."""
            ids = set()
            # There should be 118218 stars, with HIP ids in order.
            for s in stars():
                self.failIf(s[0] in ids)
                ids.add(s[0])
            self.assertEqual(len(ids), 118217)

        def testMagnitude(self):
            """Test the magnitude cutoff argument."""
            # There are only 15 stars with magnitude greater than 1.
            self.assertEqual(len(list(stars(1))), 15)

        def testDecTop(self):
            # There are only 8 stars in the top degree of celestial
            # north pole.
            self.assertEqual(len(list(stars(20,math.pi,math.radians(89)))),8)

        def testDec(self):
            # There are 939 stars with declination between 15 and 16 degrees.
            self.assertEqual(len(list(stars(
                20,math.radians(16), math.radians(15)))), 939)

        def testMagnitudeAndDec(self):
            # Moreover, there is exactly *ONE* star in the top degree
            # of the celestial north pole with magnitude at most 2,
            # Polaris.
            p = list(stars(2,math.pi,math.radians(89)))
            self.assertEqual(len(p), 1)
            self.assertEqual(p[0][0], 11767)

    # Run the tests.
    unittest.main()
