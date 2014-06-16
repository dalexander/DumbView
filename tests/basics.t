
  $ CMPH5=`python -c 'import pbcore.data as D; print D.getCmpH5()'`
  $ VARIANTS_GFF=`python -c 'import pbcore.data as D; print D.getGff3()'`
  $ LAMBDA_REF=`python -c 'import pbcore.data as D; print D.getLambdaFasta()'`
  $ DV="dumbview -N --color=off"

  $ export TERM=dumb

  $ [ -f $LAMBDA_REF   ] || echo "MISSING FILE"
  $ [ -f $VARIANTS_GFF ] || echo "MISSING FILE"
  $ [ -f $CMPH5        ] || echo "MISSING FILE"


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
  



  $ $DV $CMPH5 -w1:10 -W20
  lambda_NEB3011:1-19
                               
            ====+====|====+====
         0  GGGC-GCGACCTCGC-GGT
         1  GGGCGGCGA-CTCGCGGGT
  



  $ $DV $CMPH5 -w1:30900 -W20
  lambda_NEB3011:30890-30909
                  30900         
            |====+====|====+====
        63  GG---TGCTGCG-A-GGCGT
        64  -GC-AT-CTG-GAAGGGCGT
        65  GG--A-GCTGCGAAGGGCGT
        66  --CAATGCTGCGAA-GGCG-
        67  -GCAATGCTGC-AAGGGCGT
  



  $ $DV $CMPH5 -w1:30900 -W60
  lambda_NEB3011:30870-30929
                  30880               30900               30920         
            |====+====|====+====|====+====|====+====|====+====|====+====
        63  A-TTCAG-ATT-GCCTGACGGG---TGCTGCG-A-GGCGTTTTCCTG-TGAGGTGTCATT
        64  ATTTCAGAA-T-GCCTGAC-GGC-AT-CTG-GAAGGGCGTTTTCCTGCTGAGG-GTCATT
        65  ATTTCAGAATTAGCCTGACGGG--A-GCTGCGAAGGGCGTTTTCCT-CTGAGG-GTCATT
        66  ATTTCAGAA--AGCCTGAC--GCAATGCTGCGAA-GGCG-TTTCCTGCTGAGGTG--ATT
        67  A-TTCAGAATTAGCCTGACG-GCAATGCTGC-AAGGGCG-TTTCCTGC-GAGGTGT-ATT
  