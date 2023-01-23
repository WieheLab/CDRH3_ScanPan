#!/bin/bash

V=V3

weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f madeSeq.fasta -o madeSeq.${V}.pdf --stack-width 500
weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f Briney_CDR3.FUD.fasta -o Briney.${V}.pdf --stack-width 500
weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f igor_generated.fasta -o ior_generated.${V}.pdf --stack-width 500
weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f DH270_generated_seqs.fasta -o DH270_generated.${V}.pdf --stack-width 500
weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f regened_Logo.fasta -o generated_logo.${V}.pdf --stack-width 500
weblogo -S 1 --units probability -c chemistry -n 150 --format pdf --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f CH235_gen_seq.fasta -o CH235_generated.${V}.pdf


#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f madeSeq.fasta -o madeSeq.${V}.pdf --stack-width 500
#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f Briney_CDR3.FUD.fasta -o Briney.${V}.pdf --stack-width 500
#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f igor_generated.fasta -o ior_generated.${V}.pdf --stack-width 500
#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f DH270_generated_seqs.fasta -o DH270_generated.${V}.pdf --stack-width 500
#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f regened_Logo.fasta -o generated_logo.${V}.pdf --stack-width 500
#weblogo -S 1 --units probability -c chemistry -n 150 --format pdf -P "" --color red EDQNKRST 'neg'  --color blue HFWY '' --color black ACMPILVG '' -a ACDEFGHIKLMNPQRSTVWXY -f CH235_gen_seq.fasta -o CH235_generated.${V}.pdf

#Polar: EDQNKRST in red
#Aromatic: HFWY in blue
#Non-Polar: ACMPILVG in black

#--color red ED 'neg' --color purple QN '' --color blue KHR '' --color green SGCYT '' --color black VLPFWAIM ''
