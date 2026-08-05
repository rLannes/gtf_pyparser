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
        self.assertFalse(inter.overlaps({"start": 200, "end": 220, "strand": "+"}))
        self.assertFalse(inter.overlaps({"start": 80, "end": 100, "strand": "+"}))

        self.assertTrue(inter.overlaps({"start": 200, "end": 220, "strand": "+"}, closed=True))
        self.assertTrue(inter.overlaps({"start": 80, "end": 100, "strand": "+"}, closed=True))

    def test_classify_position(self):
        # gene span [100, 200), one exon [120, 150) inside it
        interval = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
        exon = gtf_pyparser.Interval("8", 120, 150, "+", ".", {})
        tr = gtf_pyparser.Transcript(
            transcript_id="t1", transcript_symbol=None, interval=interval,
            features={"exon": [exon]},
        )

        self.assertEqual(tr.classify_position(100, "+"), "geneStart")
        self.assertEqual(tr.classify_position(200, "+"), "geneEnd")
        self.assertIsNone(tr.classify_position(250, "+"))
        self.assertIsNone(tr.classify_position(99, "+"))
        self.assertEqual(tr.classify_position(130, "+"), "exon")
        self.assertEqual(tr.classify_position(120, "+"), "junctionAcceptor")
        self.assertEqual(tr.classify_position(150, "+"), "junctionDonnor")
        self.assertEqual(tr.classify_position(160, "+"), "intron")

    def test_intersect_strand(self):
        inter = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
        self.assertFalse(inter.overlaps({"start": 180, "end": 200, "strand": "-"}))
        self.assertTrue(inter.overlaps({"start": 180, "end": 200, "strand": "-"}, strand_aware=False))

    def test_contains(self):
        inter = gtf_pyparser.Interval("8", 100, 200, "+", ".", {})
        self.assertFalse(inter.contains(position=80, strand="+"))
        self.assertFalse(inter.contains(position=220, strand="+"))
        self.assertTrue(inter.contains(position=100, strand="+"))
        self.assertTrue(inter.contains(position=150, strand="+"))
        self.assertFalse(inter.contains(position=200, strand="+"))

        self.assertTrue(inter.contains(position=200, strand="+", closed=True))
        self.assertFalse(inter.contains(position=80, strand="+", closed=True))


if __name__ == '__main__':
    unittest.main()
