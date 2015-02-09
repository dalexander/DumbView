
  $ CMPH5=`python -c 'import pbcore.data as D; print D.getCmpH5()'`
  $ VARIANTS_GFF=`python -c 'import pbcore.data as D; print D.getGff3()'`
  $ LAMBDA_REF=`python -c 'import pbcore.data as D; print D.getLambdaFasta()'`
  $ DV="dumbview -N --color=off"

  $ export TERM=dumb

  $ [ -f $LAMBDA_REF   ] || echo "MISSING FILE"
  $ [ -f $VARIANTS_GFF ] || echo "MISSING FILE"
  $ [ -f $CMPH5        ] || echo "MISSING FILE"


First test, specifying location of cmp.h5, ref

  $ $DV -r $LAMBDA_REF $CMPH5 $VARIANTS_GFF
  deletion lambda_NEB3011 30890 30890 (G > .) 25
                           30900
       Ref  AGCCTGACGGGCAATGCTGC
            ====+====|====+====|
        63  -GCCTGACGGG---TGCTGC
        64  -GCCTGAC-GGC-AT-CTG-
        65  AGCCTGACGGG--A-GCTGC
        66  AGCCTGAC--GCAATGCTGC
        67  AGCCTGACG-GCAATGCTGC
  
  insertion lambda_NEB3011 30924 30924 (. > G) 25
             30920              
       Ref  TGCTGAGGTGTCATTGAACA
            +====|====+====|====
        63  TtGTGgAGGTGgTCATTGAtACA
        64  TGCTGcAGGGTCATTGAACA
        65  TCcTGAGGGgTaCATTGAAccC
        66  TGCTGAaGGTGATTGAACA
        67  TGgCGAGGTGTATTAACA
  
Now let it root through the GFF to find the supposed location of the
ref and alignment files.  We have to construct the GFF to point to the
pbcore cmp.h5, reference in their deployed location.

  $ cat > variants.gff.in <<EOF
  > ##gff-version 3
  > ##pacbio-variant-version 2.1
  > ##date Sat Mar 22 12:16:13 2014
  > ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  > ##source GenomicConsensus 0.8.0
  > ##source-commandline /Users/dalexander/.virtualenvs/VE/bin/variantCaller.py --algorithm=plurality -q20 -x5 pbcore/data/aligned_reads_1.cmp.h5 -r /Users/dalexander/Data/lambdaNEB.fa -o /tmp/v.gff
  > ##source-alignment-file ALIGNMENTS
  > ##source-reference-file REFERENCE
  > ##sequence-region lambda_NEB3011 1 48502
  > lambda_NEB3011	.	insertion	30924	30924	.	.	.	reference=.;variantSeq=G;frequency=2;coverage=5;confidence=25
  > EOF

  $ sed -e "s|ALIGNMENTS|$CMPH5|g" -e "s|REFERENCE|$LAMBDA_REF|g"  variants.gff.in > variants.gff


  $ $DV variants.gff
  insertion lambda_NEB3011 30924 30924 (. > G) 25
             30920              
       Ref  TGCTGAGGTGTCATTGAACA
            +====|====+====|====
        63  TtGTGgAGGTGgTCATTGAtACA
        64  TGCTGcAGGGTCATTGAACA
        65  TCcTGAGGGgTaCATTGAAccC
        66  TGCTGAaGGTGATTGAACA
        67  TGgCGAGGTGTATTAACA
  
