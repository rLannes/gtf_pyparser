import gtf_pyparser


import unittest

class TestStringMethods(unittest.TestCase):

    def test_overlaps(self):
        inter = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
        self.assertTrue(inter.overlaps({"start": 80, "end": 120, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 100, "end": 120, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 100, "end": 200, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 120, "end": 150, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 180, "end": 220, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 180, "end": 200, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 200, "end": 220, "strand": "+"}))
        self.assertTrue(inter.overlaps({"start": 80, "end": 100, "strand": "+"}))

        self.assertFalse(inter.overlaps({"start": 200, "end": 220, "strand": "+"}, semi_closed=True))
        self.assertFalse(inter.overlaps({"start": 80, "end": 100, "strand": "+"}, semi_closed=True))

       
def test_intersect_strand(self):   
        inter = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
        self.assertFalse(inter.overlaps({"start": 180, "end": 200, "strand": "-"}))
        self.assertTrue(inter.overlaps({"start": 180, "end": 200, "strand": "-"}, strand_aware=False))

def test_contains(self):
    inter = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
    self.assertFalse(self.assertTrue(inter.contains(pos=80, strand="+")))
    self.assertFalse(self.assertTrue(inter.contains(pos=220, strand="+")))
    self.assertTrue(self.assertTrue(inter.contains(pos=200, strand="+")))
    self.assertTrue(self.assertTrue(inter.contains(pos=100, strand="+")))
    self.assertTrue(self.assertTrue(inter.contains(pos=50, strand="+")))
    self.assertFalse(self.assertTrue(inter.contains(pos=200, strand="+")), semi_closed=True)
    self.assertTrue(self.assertTrue(inter.contains(pos=100, strand="+")), semi_closed=True)



if __name__ == '__main__':
    unittest.main()
