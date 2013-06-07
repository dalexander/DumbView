### *DumbView:* a PacBio genome browser for command-line nerds

![Screenshot](!https://raw.github.com/dalexander/DumbView/master/screenshot.png)

Requirements: ``pbcore``, ``GenomicConsensus``

Install: ``python setup.py install``

Use:

  - View alignments to a window of the reference:
    ``dumbview aligned_reads.cmp.h5 [-r reference.fasta] -w<contigName>:<start>-<end>``
  - View alignments around all variant calls in a ``variants.gff`` file:
    ``dumbview variants.gff``
  - To enable "unaligned" mode (useful for eyeballing insertion variants) use ``-u``

Notes:
  - DumbView expects a FASTA index (.fai) file in the same directory
    as the reference FASTA.  Make it with
    ``samtools faidx reference.fasta``.
  - Coordinates are 1-based, in agreement with the "GFF" convention
    and most other genome browsers.
  - The consensus track shown at the bottom presently uses the Quiver
    "NoQVs" model, so it may not be the same consensus as you get from
    running Quiver.  This will be fixed soon.


This is an unsupported tool and the code is horrendous.  But I find it useful.

David Alexander, 2013
