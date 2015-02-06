db=Human-Uniprot-Swissprot.fa
query=NGLY1.fa
out='NGLY1-Human-Uniprot-Swissprot-blat'
blat $db $query  -t=prot -q=prot -out=blast8 $out
