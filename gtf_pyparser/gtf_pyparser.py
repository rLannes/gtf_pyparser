import copy
from dataclasses import dataclass, field, asdict
import re
import logging
import sys


log = logging.getLogger(__name__)
# usage logging.getLogger("gtf").setLevel(logging.WARNING)      # silence gtf

@dataclass(frozen=True)
class Interval:
    """
    Immutable genomic coordinate record.

    Represents a single contiguous region on a chromosome. Coordinates follow
    0-based half-open convention (start inclusive, end exclusive), consistent
    with BED/BAM format and the conversion applied during GTF parsing.

    Attributes
    ----------
    chr_ : str
        Chromosome or sequence name (e.g. ``"chr1"``). chr_ because chr is a reserved python keyword
    start : int
        0-based start coordinate (inclusive).
    end : int
        0-based end coordinate (exclusive).
    strand : str
        Strand of the feature, typically ``"+"`` or ``"-"``. May be ``"."``
        for unstranded features.
    phase : int or str
        Reading frame phase (0, 1, or 2). Set to ``"."`` for features where
        phase is not applicable, such as introns.
    attribute : dict
        Key-value pairs parsed from the GTF attribute column, or a small
        derived dict for synthetic intervals such as introns.
    """
    chr_: str 
    start: int
    end: int
    strand: str
    phase: int = "."
    attribute: dict[str, str] = field(default_factory=dict)


    def __repr__(self):
        return "{}:{}-{}({})".format(self.chr, self.start, self.end, self.strand)

    def from_dict(dict):
        """
        chr, start, end, strand must be present or it will return a KeyError
        will also look for phase and attribute
        """
        chr_ = dict["chr"]
        start = int(dict["start"])
        end = int(dict["end"])
        strand = dict["strand"]
        phase = dict.get("phase", ".")
        attr =  dict.get("attribute")

        if attr:
            return Interval(chr_, start, end, strand, phase, attr)
        else:
            return Interval(chr_, start, end, strand, phase)
    
    @property
    def to_dict(self):
        return asdict(self)

    @property
    def length(self):
        return self.end - self.start

    @property
    def position(self):
        """
        return a dict that only have genomic position information
        """
        return {"chr": self.chr, "start": self.start, "end": self.end, "strand": self.strand}
    

    def overlaps(self, other, strand_aware=True, semi_closed=False):
        if strand_aware and self.strand != other["strand"]:
            return False
        if not semi_closed:
            if (self.end < other["start"]) or (self.start > other["end"]):
                return False
        else:
            if (self.end <= other["start"]) or (self.start >= other["end"]):
                return False
        return True
    
    def contains(self, position, strand, strand_aware=True, semi_closed=False):
        if strand_aware and self.strand != strand:
            return False
        if semi_closed:
            if self.start <= position < self.end:
                return True
        else:
            if self.start <= position <= self.end:
                return True
        return False
    
    def eq_pos(self, other):
        return self.chr == other.chr and self.start == other.start and self.end == other.end and self.strand == other.strand
    
    @property
    def chr(self):
        return self.chr_

@dataclass
class Transcript:
    """
    A single transcript isoform, with its constituent genomic features.

    Groups all GTF records sharing a common transcript_id under one object.
    The transcript's genomic span (start, end) is maintained incrementally
    as features are added via Gene._parse_transcript_line, and always
    reflects the union of all feature intervals.

    Attributes
    ----------
    transcript_id : str
        Unique identifier for the transcript (e.g. an Ensembl transcript ID).
    transcript_symbol : str or None
        Human-readable symbol for the transcript, if available.
    interval : Interval
        Current genomic span of the transcript. Replaced (not mutated) when
        new features extend the known boundaries.
    features : dict[str, list[Interval]]
        Mapping of feature type (e.g. ``"exon"``, ``"CDS"``) to a list of
        corresponding Interval records in the order they were encountered.
    """
    transcript_id: str
    transcript_symbol: str
    interval: Interval
    features: dict[str, list[Interval]] = field(default_factory=dict)
    
    # TODO 
    # def to_dict():
    # def from_dict():
    @property
    def exons(self):
        return self.features.get("exon")
    

    def __repr__(self):
        to_p = "Transcript: transcript_id {}; transcript_symbol: {}, length: {}\n".format(self.transcript_id, self.transcript_symbol, self.interval.length)
        to_p += "{}\n".format(self.interval)

        if (exon := self.exons_positions):
            to_p += "number_exon: {}\n".format(len(exon))
        return to_p

    @property
    def cds(self): 
        if (cds := self.features.get("CDS")):
            return cds 
        if (cds := self.features.get("cds")):
            return cds 
    
    @property
    def mrna(self):
        if (mrna := self.features.get("mRNA")):
            return mrna 
        if (mrna := self.features.get("transcript")):
            return mrna 

    @property
    def length(self):
        return self.interval.length
    
    @property
    def utr_5p(self):
        if (utr := self.features.get("5UTR")):
            return utr # flybase
        if (utr := self.features.get("five_prime_utr")): # ensembl
            return utr

    @property
    def utr_3p(self):
        if (utr := self.features.get("3UTR")  ):
            return utr # flybase
        if (utr := self.features.get("three_prime_utr")): # ensembl
            return utr
        
    @property
    def exons_positions(self):
        if "exon" in self.features:
            return [x.position for x in self.features.get("exon")]
        else: 
            return None

    @property
    def attribute(self):
        return self.interval.attribute

    @property
    def start(self):
        return self.interval.start
    
    @property
    def chr(self):
        return self.interval.chr
    
    @property
    def end(self):
        return self.interval.end

    @property
    def phase(self):
        return self.interval.phase
        
    @property
    def intron(self):
        if len(self.exons) > 1:
            return get_intron(self.exons)
        else:
            return []
        
    @phase.setter
    def phase(self, value):
        self.interval.phase = value

                
    def classify_position(self, position, strand, strand_aware=True, semi_closed=False):
        # -> "None", "exon", "intron", "junctionDonnor", "junctionAcceptor", "geneStart", "geneEnd"
        if not self.interval.contains(position=position, strand=strand, strand_aware=strand_aware, semi_closed=semi_closed):
            return None
        if self.interval.start == position:
            return "geneStart"
        if  self.interval.end == position:
            return "geneEnd"

        # do they interesect with exon:
        for exon in self.exons:
            if exon.contains(position=position, strand=strand, strand_aware=strand_aware, semi_closed=semi_closed):
                if exon.start == position:
                    if "strand" == "+":
                        return "junctionAcceptor"
                    return "junctionDonnor"
                elif exon.end == position:
                    if "strand" != "+":
                        return "junctionAcceptor"
                    return "junctionDonnor"
                else:
                    return "exon"
        return "intron"



@dataclass
class Gene:
    """
    A gene locus, containing one or more transcript isoforms.

    Wraps a genomic Interval with a collection of Transcript objects and
    convenience properties for commonly accessed attributes. The gene span
    is set directly from ``gene``-type GTF records; transcript spans are
    managed independently by _parse_transcript_line.

    Attributes
    ----------
    gene_id : str
        Primary identifier for the gene (typically an Ensembl gene ID).
    symbol : str
        Human-readable gene symbol, resolved from ``gene_symbol``,
        ``gene_name``, or gene_id in that order.
    interval : Interval
        Genomic span of the gene as recorded in the GTF gene record.
    transcripts : dict[str, Transcript]
        Mapping of transcript_id to Transcript objects associated with
        this gene.
    """
    gene_id: str
    symbol: str
    interval: Interval
    transcripts: dict[str, Transcript] = field(default_factory=dict)
    
    # TODO 
    # def to_dict():
    # def from_dict():

    @property
    def length(self):
        return self.interval.length

    @property
    def attribute(self):
        return self.interval.attribute
    
    @property
    def transcript_names(self):
        return list(self.transcripts.keys())

    @property
    def transcript_length(self):
        return len(self.transcripts)

    @property
    def biotype(self):
        return self.interval.attribute.get("biotype", None)
    
    @property
    def has_transcript(self):
        return len(self.transcripts) > 0
    
    @property
    def has_exon(self):
        if not self.has_transcript:
            return False
        for tr_id, tr_item in self.transcripts.items():
            if tr_item.exons:
                return True
        return False

    def __repr__(self):
        to_p = "Gene: gene_id: {}; gene_symbol: {}\n".format(self.gene_id, self.symbol)
        to_p += "{}\n".format(self.interval)
        for tr_id, tr in self.transcripts.items():
            to_p += "{}\n".format(tr)

        return to_p

    
    @property
    def start(self):
        return self.interval.start
    
    @property
    def end(self):
        return self.interval.end

    @property
    def chr(self):
        return self.interval.chr

    @property
    def phase(self):
        return self.interval.phase

    def classify_position(self, position, strand, strand_aware=True, semi_closed=False):
        # look at all transcript and return all hit
        # return -> dict[key] -> "None", "exon", "intron", "junctionDonnor", "junctionAcceptor", "geneStart", "geneEnd"
        res = {}
        for transcript_key, transcript_item in self.transcripts.items():
            res[transcript_key] =  transcript_item.classify_position(position=position, strand=strand, strand_aware=strand_aware, semi_closed=semi_closed)
        return res

    # private helper method that delegate parsing transript out of the main loop
    def _parse_transcript_line(self, transcript_id, transcript_symbol,  interval, feature):
        """
        Register a single GTF feature record under the appropriate transcript.

        If the transcript does not yet exist, it is initialised with a deep copy
        of the provided interval as its initial span. The feature interval is
        appended to the transcript's feature list for the given feature type.
        The transcript's genomic span is then updated if the new interval extends
        beyond the current known boundaries, replacing the interval object rather
        than mutating it in place.

        Parameters
        ----------
        transcript_id : str
            Identifier of the transcript this record belongs to.
        transcript_symbol : str or None
            Human-readable symbol for the transcript, passed through on
            first creation only.
        interval : Interval
            The genomic interval of this feature record, as parsed from the GTF.
        feature : str
            Feature type string from GTF column 3 (e.g. ``"exon"``, ``"CDS"``,
            ``"start_codon"``).

        Notes
        -----
        This is an internal helper called exclusively from gtf_to_dict and is not
        part of the public API.
        """
        if transcript_id not in self.transcripts:
            self.transcripts[transcript_id] = Transcript(transcript_id=transcript_id, transcript_symbol=transcript_symbol,
                                                       interval=copy.deepcopy(interval))
        
        if feature not in self.transcripts[transcript_id].features:
            self.transcripts[transcript_id].features[feature] = []
        self.transcripts[transcript_id].features[feature].append(copy.deepcopy(interval))
        tr_i = self.transcripts[transcript_id].interval

        if interval.start < self.transcripts[transcript_id].start:
            interval = Interval(tr_i.chr, start=interval.start, end=tr_i.end, strand=tr_i.strand, phase=".", attribute=tr_i.attribute)
            self.transcripts[transcript_id].interval = interval
        if interval.end > self.transcripts[transcript_id].end:
            interval = Interval(tr_i.chr, start=tr_i.start, end=interval.end, strand=tr_i.strand, phase=".", attribute=tr_i.attribute)
            self.transcripts[transcript_id].interval = interval


REG = re.compile(r'(\w+)\s+"([^"]*)"')
def get_attr(string, reg=REG):
    """
    Parse a GTF attribute string into a key-value dictionary.

    GTF attribute fields consist of semicolon-separated tag-value pairs, where each
    pair has the form ``tag "value"`` or ``tag value``. Surrounding whitespace and
    double quotes are stripped from values.

    Parameters
    ----------
    string : str
        The raw attribute field from column 9 of a GTF record, e.g.
        ``'gene_id "ENSG00000001"; transcript_id "ENST00000001";'``.

    Returns
    -------
    dict[str, str]
        Mapping of attribute tag to value. Empty tokens produced by trailing
        semicolons are silently skipped.

    Raises
    ------
    Exception
        Re-raises any parsing exception after logging the offending token list
        at ERROR level.
    """
    dico = {}

    for (key, value) in reg.findall(string):
        try:
            this_key = key.strip()
            if this_key in dico:
                v = dico[this_key]
                if not isinstance(v, list):
                    dico[this_key] = [v]
                dico[this_key].append(value)
            else:
                dico[this_key] = value.replace('"', "").strip() 
        except:
            logging.error("failed to parse line {}".format(string))
            raise
    return dico


class Gtf():
    
    def __init__(self):
        self.dict = {}

    def __setitem__(self, key, value):
        self.dict[key] = value
    
    def __getitem__(self, key):
        return self.dict[key]
    
    def __delitem__(self, key):
        del self.dict[key]

    def __contains__(self, key):
        return key in self.dict

    def __len__(self):
        return len(self.dict )
    
    def __iter__(self):
        return (key for key in self.dict)
    
    def keys(self):
        return self.dict.keys()

    def values(self):
        return self.dict.values()

    def items(self):
        return self.dict.items()

    # TODO 
    # def to_dict():
    # def from_dict():
    # def add_intervalTree add an intervall tree alowing to recover ggene_id by position? or Gen obj # required import interval Tree import intervaltree as it
    # becasue the lib is no dependancy the import will be maden only if user use this function.
    # number genes, avergae gene length, gene by chromosome etc....


# Should main dico be a dataclass to or an object I could implement getitem setitem on it so it behave like a dict but I will have method on top of it

def gtf_to_dict(gtf_file, primary_key = "gene_id"):
    """
    Parse a GTF file into a dictionary of Gene objects keyed by a primary attribute.

    Reads the file line by line, skipping blank lines and comment lines (starting
    with ``#``). For each record, coordinates are converted from 1-based inclusive
    (GTF convention) to 0-based half-open. Gene records initialise a new Gene entry
    if one does not already exist. Transcript and feature records are delegated to
    Gene._parse_transcript_line, which maintains per-transcript interval bounds
    incrementally.

    Parameters
    ----------
    gtf_file : str or path-like
        Path to the GTF file to parse.
    primary_key : str, optional
        Attribute tag used to group records into genes. Defaults to ``"gene_id"``.
        The tag must be present in every record's attribute field.

    Returns
    -------
    dict[str, Gene]
        Mapping of primary key value to Gene object. Each Gene contains all
        transcripts and features encountered for that key in the file.

    Raises
    ------
    KeyError
        If a record's attribute field does not contain the primary key.
    Exception
        Re-raises any line-level parsing error after logging the offending line
        at ERROR level.

    Notes
    -----
    Gene symbol is resolved by looking up ``gene_symbol`` first, then falling back
    to ``gene_name``, then to the primary key value itself.
    Only records carrying a ``transcript_id`` attribute are added to transcript
    feature lists; bare gene-level records are used solely to set gene coordinates.
    Additionaly each Interval carries "feature_": <feature type / 3rd field of the file> attribute.
    """
    dico = Gtf()

    log.info("reading gtf file {} into a dict".format(gtf_file))
    with open(gtf_file) as f_in:
        """atm reads only genes"""
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            spt = line.strip().split("\t")
            try:
                chr_ = spt[0]
                start = int(spt[3]) - 1 # gtf are 1 based convert to 0 based
                end = int(spt[4])
                strand =  spt[6]            
                phase = spt[7]
                attr = get_attr(spt[-1])
                type_ = spt[2]
                attr["feature_"] = type_
                try:
                    gene_id = attr[primary_key]
                except:
                    log.error("failed to recover gene_id {} {}".format(attr, spt))
                    raise
                gene_symbol = attr.get("gene_symbol", gene_id)
                if "gene_symbol" not in attr:
                    gene_symbol = attr.get("gene_name", gene_id)
                
            except:
                log.error("failed to parse line: {} {} {}".format(attr, line, spt))
                raise
            
            this_interval = Interval(chr_=chr_, start=start, end=end, strand=strand, phase=phase, attribute=attr)

            if gene_id not in dico:
                
                this = Gene(gene_id=gene_id, symbol=gene_symbol,
                     interval=copy.deepcopy(this_interval))
                     
                dico[gene_id] = this
                
                
            if type_ == "gene":
                dico[gene_id].interval = copy.deepcopy(this_interval)
            
            # if not a transcript pass
            if (transcript_id  := attr.get("transcript_id", None)):
                transcript_symbol = attr.get("transcript_symbol", None)
                if not transcript_symbol:
                    transcript_symbol = attr.get("transcript_name", None)
                dico[gene_id]._parse_transcript_line(transcript_id, transcript_symbol, this_interval, type_)


    return dico



def get_intron(exons: list[Interval]):
    """
    Derive intron intervals from a list of exon intervals belonging to a single transcript.

    Introns are inferred as the gaps between consecutive exons after sorting by start
    position. Each returned Interval spans from the end of one exon to the start of the
    next, and carries an ``intron_n`` attribute reflecting biological ordering along the
    transcript: introns are numbered starting from 1 on the plus strand, and from
    n_introns down to 1 on the minus strand. Additionaly each Interval carries "feature_": "intron" attribute

    Parameters
    ----------
    exons : list[Interval]
        Exon intervals for a single transcript. All intervals are assumed to share the
        same chromosome and strand. Coordinates must satisfy start < end (0-based,
        half-open). The list need not be pre-sorted.

    Returns
    -------
    list[Interval]
        Intron intervals in genomic order (sorted by start position), each with:
        - ``phase`` set to ``"."``
        - ``attribute`` containing a single key ``intron_n`` (int)
        Returns an empty list if fewer than two exons are provided.

    Notes
    -----
    The strand of each returned intron is taken from the first exon after sorting, which
    is assumed to be consistent across all input exons.
    
    
    """
    if len(exons) <= 1: 
        return []
    
    exon_sorted = sorted(exons, key=lambda x: x.start)
    n_introns = len(exon_sorted) - 1

    this_n = 1 if exon_sorted[0].strand == "+" else n_introns
    strand = exon_sorted[0].strand
    introns = []

    for i in range(len(exon_sorted) - 1 ):
        e = exon_sorted[i]
        e0 = exon_sorted[i].end
        e1 = exon_sorted[i+1].start

        attr = {}
        attr["intron_n"] =  this_n
        attr["feature_"] = "intron"
        

        this_n += 1 if exon_sorted[0].strand == "+" else -1
        
        introns.append(Interval(chr_=e.chr, start=e0, end=e1, strand=strand, phase=".", attribute=attr))

    return introns

