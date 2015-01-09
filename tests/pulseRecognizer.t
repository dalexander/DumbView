
  $ CMPH5=`python -c 'import pbcore.data as D; print D.getCmpH5()'`
  $ VARIANTS_GFF=`python -c 'import pbcore.data as D; print D.getGff3()'`
  $ LAMBDA_REF=`python -c 'import pbcore.data as D; print D.getLambdaFasta()'`
  $ DV="dumbview -N --color=off"

  $ export TERM=dumb

  $ [ -f $LAMBDA_REF   ] || echo "MISSING FILE"
  $ [ -f $VARIANTS_GFF ] || echo "MISSING FILE"
  $ [ -f $CMPH5        ] || echo "MISSING FILE"

  $ $DV --pulseRecognizer -r $LAMBDA_REF $CMPH5 $VARIANTS_GFF
  deletion lambda_NEB3011 30890 30890 (G > .) 25
                           30900
       Ref  AGCCTGACGGGCAATGCTGC
            ====+====|====+====|
        63  -GCCTGACGGG---TGCTGC
        64  -GCCTGAC-GGC-AT-CTG-
        65  AGCCTGACGGG--A-GCTGC
        66  AGCCTGAC--GCAATGCTGC
        67  AGCCTGACG-GCAATGCTGC
  
    Link for pulse recognizer: <a href=./variant-0.csv>variant-0.csv</a>
  
  
  insertion lambda_NEB3011 30924 30924 (. > G) 25
             30920              
       Ref  TGCTGAGGTGTCATTGAACA
            +====|====+====|====
        63  TtGTGgAGGTGgTCATTGAtACA
        64  TGCTGcAGGGTCATTGAACA
        65  TCcTGAGGGgTaCATTGAAccC
        66  TGCTGAaGGTGATTGAACA
        67  TGgCGAGGTGTATTAACA
  
    Link for pulse recognizer: <a href=./variant-1.csv>variant-1.csv</a>
  
  





  $ cat variant-0.csv
  Movie,HoleNumber,rStart,rEnd
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2953,2963
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2350,2362
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,1572,1581
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2154,2163
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,3139,3152

  $ cat variant-1.csv
  Movie,HoleNumber,rStart,rEnd
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2915,2928
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2390,2401
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,1608,1620
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,2122,2132
  m110818_075520_42141_c100129202555500000315043109121112_s1_p0,2008,3176,3186
