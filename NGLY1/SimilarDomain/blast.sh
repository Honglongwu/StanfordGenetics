db=Human-Uniprot-Swissprot.fa
query=PNG1.fa
out='PNG1-Uniprot-Swissprot-blast'
blastp  -db $db  -query $query -out $out -outfmt 6  
#blastp  -db $db  -query $query -out $out 
