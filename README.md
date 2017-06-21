# imgt2aa
extract aligned amino acid sequences from IMGT/HLA


# Prerequsites

Need to first pull the hla.xml file

```
$ curl ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip >hla.xml
.zip
$ unzip hla.xml.zip
```

# Run

perl imgt2aa.pl >DPB1.db


# How it works

The hla.xml file contains nucleotide sequences and cDNA coordinates for HLA alleles.  The nucleotides for Exon2 are extracted and parsed.  Then the cDNA coordinates are used to offset the AA sequence such that the first postion in the string (possibly "*") corresponds to AA 1 in the mature protine


# TODO

* Complete implemtnation for all loci. Only implemented for DPB1
* Incoporate other exons.  Only implemented for exon 2



Martin Maiers
<mmaiers@nmdp.org>
