import gtf_pyparser


import unittest

try:
    import intervaltree
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False

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
        # boundary check matches the raw (exclusive) end value, 200 -- not the last real base (199)
        self.assertEqual(tr.classify_position(200, "+"), "geneEnd")
        self.assertEqual(tr.classify_position(199, "+"), "intron")  # last real base, not itself the boundary
        self.assertIsNone(tr.classify_position(201, "+"))  # past the closed [start, end] span
        self.assertIsNone(tr.classify_position(250, "+"))
        self.assertIsNone(tr.classify_position(99, "+"))
        self.assertEqual(tr.classify_position(130, "+"), "exon")
        self.assertEqual(tr.classify_position(120, "+"), "exonAcceptor")
        # boundary check matches the raw (exclusive) end value, 150
        self.assertEqual(tr.classify_position(150, "+"), "exonDonor")
        self.assertEqual(tr.classify_position(149, "+"), "exon")  # exon's last real base, not itself the boundary
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


@unittest.skipUnless(HAS_INTERVALTREE, "requires optional intervaltree dependency")
class TestGetGenesAtPosition(unittest.TestCase):

    def _make_gtf(self):
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gtf.dict = {
            "G1": gtf_pyparser.Gene(gene_id="G1", symbol="G1",
                                     interval=gtf_pyparser.Interval("1", 100, 200, "+")),
            "G2": gtf_pyparser.Gene(gene_id="G2", symbol="G2",
                                     interval=gtf_pyparser.Interval("1", 300, 400, "-")),
            "G3": gtf_pyparser.Gene(gene_id="G3", symbol="G3",
                                     interval=gtf_pyparser.Interval("2", 100, 200, "+")),
        }
        return gtf

    def test_overlap_returns_gene_id(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 150, "end": 160})
        self.assertEqual(res, ["G1"])

    def test_no_overlap_returns_none(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 220, "end": 230})
        self.assertIsNone(res)

    def test_unknown_chromosome_returns_none(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "X", "start": 10, "end": 20})
        self.assertIsNone(res)

    def test_builds_interval_tree_lazily(self):
        gtf = self._make_gtf()
        self.assertIsNone(gtf.it)
        gtf.get_genes_at_position({"chr": "1", "start": 150, "end": 160})
        self.assertIsNotNone(gtf.it)

    def test_flanking_extends_query(self):
        gtf = self._make_gtf()
        # position just past the gene's end: no overlap without flanking...
        self.assertIsNone(gtf.get_genes_at_position({"chr": "1", "start": 205, "end": 210}))
        # ...but flanking pulls it into range.
        res = gtf.get_genes_at_position({"chr": "1", "start": 205, "end": 210}, flanking=10)
        self.assertEqual(res, ["G1"])

    def test_strand_aware_filters_by_strand(self):
        gtf = self._make_gtf()
        position = {"chr": "1", "start": 350, "end": 360}
        res_plus = gtf.get_genes_at_position({**position, "strand": "+"}, strand_aware=True)
        self.assertEqual(res_plus, [])
        res_minus = gtf.get_genes_at_position({**position, "strand": "-"}, strand_aware=True)
        self.assertEqual(res_minus, ["G2"])

    def test_multiple_overlapping_genes(self):
        gtf = self._make_gtf()
        gtf.dict["G4"] = gtf_pyparser.Gene(gene_id="G4", symbol="G4",
                                            interval=gtf_pyparser.Interval("1", 150, 250, "+"))
        res = gtf.get_genes_at_position({"chr": "1", "start": 160, "end": 170})
        self.assertEqual(sorted(res), ["G1", "G4"])

    # -- half-open boundary edge cases: G1 spans [100, 200) --

    def test_boundary_touching_start_not_overlapping(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 90, "end": 100})
        self.assertIsNone(res)

    def test_boundary_at_gene_start_overlaps(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 100, "end": 101})
        self.assertEqual(res, ["G1"])

    def test_boundary_touching_end_not_overlapping(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 200, "end": 210})
        self.assertIsNone(res)

    def test_boundary_last_base_overlaps(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 199, "end": 200})
        self.assertEqual(res, ["G1"])

    def test_empty_gtf_returns_none(self):
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gtf.dict = {}
        res = gtf.get_genes_at_position({"chr": "1", "start": 10, "end": 20})
        self.assertIsNone(res)

    def test_index_is_cached_not_rebuilt(self):
        gtf = self._make_gtf()
        gtf.build_intervals()
        # Added after the index was built: get_genes_at_position must reuse
        # the cached tree (self.it is already truthy) rather than rebuild it,
        # so this gene stays invisible until build_intervals runs again.
        gtf.dict["G5"] = gtf_pyparser.Gene(gene_id="G5", symbol="G5",
                                            interval=gtf_pyparser.Interval("1", 500, 600, "+"))
        res = gtf.get_genes_at_position({"chr": "1", "start": 550, "end": 560})
        self.assertIsNone(res)

    def test_strand_aware_defaults_to_plus_when_missing(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "1", "start": 150, "end": 160}, strand_aware=True)
        self.assertEqual(res, ["G1"])

        res_none = gtf.get_genes_at_position({"chr": "1", "start": 350, "end": 360}, strand_aware=True)
        self.assertEqual(res_none, [])

    def test_nested_genes_both_returned(self):
        gtf = self._make_gtf()
        gtf.dict["G_outer"] = gtf_pyparser.Gene(gene_id="G_outer", symbol="G_outer",
                                                 interval=gtf_pyparser.Interval("3", 100, 1000, "+"))
        gtf.dict["G_inner"] = gtf_pyparser.Gene(gene_id="G_inner", symbol="G_inner",
                                                 interval=gtf_pyparser.Interval("3", 400, 500, "+"))
        res = gtf.get_genes_at_position({"chr": "3", "start": 450, "end": 460})
        self.assertEqual(sorted(res), ["G_inner", "G_outer"])

    def test_only_queried_chromosome_matches(self):
        gtf = self._make_gtf()
        res = gtf.get_genes_at_position({"chr": "2", "start": 150, "end": 160})
        self.assertEqual(res, ["G3"])

    def test_flanking_needed_only_on_left_side(self):
        gtf = self._make_gtf()
        self.assertIsNone(gtf.get_genes_at_position({"chr": "1", "start": 90, "end": 95}))
        res = gtf.get_genes_at_position({"chr": "1", "start": 90, "end": 95}, flanking=10)
        self.assertEqual(res, ["G1"])


@unittest.skipUnless(HAS_INTERVALTREE, "requires optional intervaltree dependency")
class TestGtfClassifyPosition(unittest.TestCase):

    def _make_gtf(self):
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gtf.dict = {
            "G1": gtf_pyparser.Gene(gene_id="G1", symbol="G1",
                                     interval=gtf_pyparser.Interval("1", 100, 200, "+")),
        }
        return gtf

    def test_returns_transcript_classification_for_single_gene(self):
        gtf = self._make_gtf()
        gtf.dict["G1"].transcripts["G1.t1"] = gtf_pyparser.Transcript(
            transcript_id="G1.t1", transcript_symbol=None,
            interval=gtf_pyparser.Interval("1", 100, 200, "+"),
            features={"exon": [gtf_pyparser.Interval("1", 140, 160, "+")]},
        )

        result = gtf.classify_position({"chr": "1", "start": 150, "strand": "+"}, strand_aware=False)

        self.assertEqual(len(result), 1)
        gene_id, classification = result[0]
        self.assertEqual(gene_id, "G1")
        self.assertEqual(classification, {"G1.t1": "exon"})

    def _gene_with_transcripts(self, gene_id, start, end, strand, transcripts):
        gene = gtf_pyparser.Gene(gene_id=gene_id, symbol=gene_id,
                                  interval=gtf_pyparser.Interval("1", start, end, strand))
        for tr_id, exon in transcripts.items():
            features = {"exon": [exon]} if exon else {}
            gene.transcripts[tr_id] = gtf_pyparser.Transcript(
                transcript_id=tr_id, transcript_symbol=None,
                interval=gtf_pyparser.Interval("1", start, end, strand),
                features=features,
            )
        return gene

    def test_returns_all_overlapping_genes_with_multiple_transcripts(self):
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gtf.dict = {
            # two overlapping genes at position 155, each with several
            # transcript isoforms, to check every overlapping gene surfaces
            # (not just the first one found in the interval tree) and that
            # each transcript is classified independently.
            "G1": self._gene_with_transcripts("G1", 100, 200, "+", {
                "G1.t1": gtf_pyparser.Interval("1", 140, 160, "+"),  # position falls in this exon
                "G1.t2": None,  # no exon: position falls in an intron
            }),
            "G2": self._gene_with_transcripts("G2", 150, 250, "+", {
                "G2.t1": None,
            }),
        }

        result = gtf.classify_position({"chr": "1", "start": 155, "strand": "+"}, strand_aware=False)

        self.assertEqual(len(result), 2)
        by_gene_id = dict(result)
        self.assertEqual(by_gene_id["G1"], {"G1.t1": "exon", "G1.t2": "intron"})
        self.assertEqual(by_gene_id["G2"], {"G2.t1": "intron"})

    def _run(self, gtf, position, strand_aware):
        result = gtf.classify_position(position, strand_aware=strand_aware)
        return [classification for _gene_id, classification in result]

    def _single_transcript_gtf(self, strand, exon=None):
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        features = {"exon": [exon]} if exon else {}
        gene = gtf_pyparser.Gene(gene_id="G1", symbol="G1",
                                  interval=gtf_pyparser.Interval("1", 100, 200, strand))
        gene.transcripts["G1.t1"] = gtf_pyparser.Transcript(
            transcript_id="G1.t1", transcript_symbol=None,
            interval=gtf_pyparser.Interval("1", 100, 200, strand),
            features=features,
        )
        gtf.dict = {"G1": gene}
        return gtf

    def _single_transcript_gtf_wider_gene(self, strand):
        # gene span [100, 205) is wider than the transcript's own [100, 200):
        # get_genes_at_position uses half-open matching, so it would never
        # find a gene whose own span ends exactly at the queried position
        # (see test_gene_end_unreachable_when_gene_span_matches_transcript
        # below). Widening the gene span here isolates classify_position's
        # own closed-boundary logic from that outer-lookup gap.
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gene = gtf_pyparser.Gene(gene_id="G1", symbol="G1",
                                  interval=gtf_pyparser.Interval("1", 100, 205, strand))
        gene.transcripts["G1.t1"] = gtf_pyparser.Transcript(
            transcript_id="G1.t1", transcript_symbol=None,
            interval=gtf_pyparser.Interval("1", 100, 200, strand),
        )
        gtf.dict = {"G1": gene}
        return gtf

    def test_gene_start_and_end_on_plus_strand(self):
        gtf = self._single_transcript_gtf_wider_gene("+")
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 100, "strand": "+"}, False),
                          [{"G1.t1": "geneStart"}])
        # transcript span is [100, 200): the boundary check matches the raw
        # (exclusive) end value, 200 -- not the last real base (199)
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 200, "strand": "+"}, False),
                          [{"G1.t1": "geneEnd"}])

    def test_gene_start_and_end_on_minus_strand(self):
        gtf = self._single_transcript_gtf_wider_gene("-")
        # on the minus strand geneStart/geneEnd flip relative to interval start/end
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 100, "strand": "-"}, False),
                          [{"G1.t1": "geneEnd"}])
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 200, "strand": "-"}, False),
                          [{"G1.t1": "geneStart"}])

    def test_gene_end_reachable_when_gene_span_matches_transcript(self):
        # gene span == transcript span [100, 200): the interval-tree lookup
        # still finds the gene (it pads the query by one base on each side),
        # and 200 -- the boundary (exclusive) end value -- is correctly
        # labeled geneEnd.
        gtf = self._single_transcript_gtf("+")
        result = gtf.classify_position({"chr": "1", "start": 200, "strand": "+"}, strand_aware=False)
        self.assertEqual(result, [("G1", {"G1.t1": "geneEnd"})])

    def test_position_past_transcript_end_is_none_even_when_gene_is_found(self):
        # position 201 is one base past the transcript's closed [100, 200]
        # span. The gene span is wider ([100, 205)), so the interval-tree
        # lookup still surfaces the gene as a candidate, but the transcript
        # itself must correctly classify 201 as outside its span rather than
        # mislabeling it geneEnd.
        gtf = self._single_transcript_gtf_wider_gene("+")
        result = gtf.classify_position({"chr": "1", "start": 201, "strand": "+"}, strand_aware=False)
        self.assertEqual(result, [("G1", {"G1.t1": None})])

    def test_exon_acceptor_and_donor_on_plus_strand(self):
        exon = gtf_pyparser.Interval("1", 150, 180, "+")
        gtf = self._single_transcript_gtf("+", exon=exon)
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 150, "strand": "+"}, False),
                          [{"G1.t1": "exonAcceptor"}])
        # 179 is the exon's real last base, not itself the boundary
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 179, "strand": "+"}, False),
                          [{"G1.t1": "exon"}])
        # the boundary check matches the raw (exclusive) end value, 180
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 180, "strand": "+"}, False),
                          [{"G1.t1": "exonDonor"}])

    def test_exon_acceptor_and_donor_on_minus_strand(self):
        exon = gtf_pyparser.Interval("1", 150, 180, "-")
        gtf = self._single_transcript_gtf("-", exon=exon)
        # on the minus strand acceptor/donor flip relative to exon start/end
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 150, "strand": "-"}, False),
                          [{"G1.t1": "exonDonor"}])
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 180, "strand": "-"}, False),
                          [{"G1.t1": "exonAcceptor"}])

    def test_none_when_transcript_does_not_reach_position(self):
        # gene span [100, 300) overlaps the query at the gene-index level, but
        # this transcript's own span only reaches [100, 200) -- it should be
        # classified None rather than "intron"/"exon".
        gtf = gtf_pyparser.Gtf(gtf_file="dummy.gtf")
        gene = gtf_pyparser.Gene(gene_id="G1", symbol="G1",
                                  interval=gtf_pyparser.Interval("1", 100, 300, "+"))
        gene.transcripts["G1.t1"] = gtf_pyparser.Transcript(
            transcript_id="G1.t1", transcript_symbol=None,
            interval=gtf_pyparser.Interval("1", 100, 200, "+"),
        )
        gtf.dict = {"G1": gene}

        self.assertEqual(self._run(gtf, {"chr": "1", "start": 250, "strand": "+"}, False),
                          [{"G1.t1": None}])

    def test_none_when_strand_mismatches_and_strand_aware(self):
        # get_genes_at_position is always called strand-unaware internally, so
        # the gene still surfaces; strand_aware=True should then make the
        # per-transcript classification None on a strand mismatch.
        gtf = self._single_transcript_gtf("+")
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 150, "strand": "-"}, True),
                          [{"G1.t1": None}])

    def test_classified_when_strand_mismatches_and_not_strand_aware(self):
        # with strand_aware=False a mismatched query strand must not suppress
        # classification -- the transcript's own strand ("+") still drives
        # geneStart/geneEnd naming, ignoring the query strand ("-").
        gtf = self._single_transcript_gtf("+")
        self.assertEqual(self._run(gtf, {"chr": "1", "start": 150, "strand": "-"}, False),
                          [{"G1.t1": "intron"}])

        gtf_wider = self._single_transcript_gtf_wider_gene("+")
        self.assertEqual(self._run(gtf_wider, {"chr": "1", "start": 100, "strand": "-"}, False),
                          [{"G1.t1": "geneStart"}])


if __name__ == '__main__':
    unittest.main()
