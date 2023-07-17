# imgt2aa
extract aligned amino acid sequences from IMGT/HLA


# Prerequisites

Need to first pull the hla.xml file

```
$ curl ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip -o hla.xml.zip
$ unzip hla.xml.zip
```

# Run

```
$ perl imgt2aa.pl 
```


Generates output files:
```
HLA-A.db
HLA-b.db
HLA-C.db
HLA-DRB1.db
HLA-DQB1.db
HLA-DPB1.db
```

Each output file is tab-delimited with:
* locus
* allele (first 2 fields)
* amino acid sequence presented such that the position in the string corresponds to the position in the mature protein



# How it works

The *hla.xml* file contains nucleotide sequences and cDNA coordinates for HLA alleles.  The nucleotides for Exon2 are extracted and parsed.  Then the cDNA coordinates are used to offset the AA sequence such that the first position in the string (possibly "\*") corresponds to AA 1 in the mature protein.


# TODO

* Complete implementation for all loci. Only implemented for A, B, C, DRB1, DQB, DPB1
* Incorporate other exons.  Only implemented for ARS exons
* Address reading frame on the 3' end ("-")
* generalize to work with GFE defined alleles (a superset of IMGT/HLA)
* move to .DAT file which has BioPython support
* Get it to work with KIR and other loci



Martin Maiers
<mmaiers@nmdp.org>
