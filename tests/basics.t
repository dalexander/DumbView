
  $ CMPH5=`python -c 'import pbcore.data as D; print D.getCmpH5()'`
  $ VARIANTS_GFF=`python -c 'import pbcore.data as D; print D.getGff3()'`
  $ DV="dumbview -N --color=off"

  $ export TERM=dumb

  $ $DV $CMPH5 $VARIANTS_GFF
  deletion lambda_NEB3011 30890 30890 (G > .) 25
                           30900
       Ref  AGCCTGACGGGCAATGCTGC
            ====+====|====+====|
        63  -GCCTGACGGG---TGCTGC
        64  -GCCTGACG-GCA-T-CTG-
        65  AGCCTGACGGG--A-GCTGC
        66  AGCCTGACG--CAATGCTGC
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
  lambda_NEB3011:0-19
                                
            |====+====|====+====
         0   GGGCG-CGACCTCGCGG-T
         1   GGGCGGCGAC-TCGCGGGT
  



  $ $DV $CMPH5 -w1:30900 -W20
  lambda_NEB3011:30890-30909
                  30900         
            |====+====|====+====
        63  GG---TGCTGCGA--GGCGT
        64  -GCA-T-CTG-GAAGGGCGT
        65  GG--A-GCTGCGAAGGGCGT
        66  --CAATGCTGCGAAGG-CG-
        67  -GCAATGCTGC-AAGGGCGT
  



  $ $DV $CMPH5 -w1:30900 -W60
  lambda_NEB3011:30870-30929
                  30880               30900               30920         
            |====+====|====+====|====+====|====+====|====+====|====+====
        63  A-TTCAG-ATT-GCCTGACGGG---TGCTGCGA--GGCGTTTTCCTG-TGAGGTGTCATT
        64  ATTTCAGAAT--GCCTGACG-GCA-T-CTG-GAAGGGCGTTTTCCTGCTGAGG-GTCATT
        65  ATTTCAGAATTAGCCTGACGGG--A-GCTGCGAAGGGCGTTTTCCT-CTGAGG-GTCATT
        66  ATTTCAGAA--AGCCTGACG--CAATGCTGCGAAGG-CG-TTTCCTGCTGAGGTG--ATT
        67  A-TTCAGAATTAGCCTGACG-GCAATGCTGC-AAGGGCGT-TTCCTGC-GAGGTGT-ATT
  
