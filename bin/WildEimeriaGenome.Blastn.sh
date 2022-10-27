

makeblastdb -in /SAN/Eimeria_Wild_Genomes/Assembly_Masurca/E64.final.genome.scf.fasta -dbtype nucl -out tmp/E164.Masurca
makeblastdb -in /SAN/Eimeria_Wild_Genomes/Assembly_Masurca/E139.final.genome.scf.fasta -dbtype nucl -out tmp/E139.Masurca



blastn -query /SAN/Susanas_den/gitProj/HMHZ/tmp/Eimeria28S.fasta -db tmp/E64.Masurca -max_target_seqs 5 -out tmp/E64_28S_blasn.txt
blastn -query /SAN/Susanas_den/gitProj/HMHZ/tmp/Eimeria28S.fasta -db tmp/E139.Masurca -max_target_seqs 5 -out tmp/E139_28S_blasn.txt

makeblastdb -in /SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa -dbtype nucl -out tmp/SilvaLSU
blastn -query /SAN/Susanas_den/gitProj/HMHZ/tmp/Eimeria28S.fasta -db tmp/SilvaLSU -max_target_seqs 5 -out tmp/LSU_28S_blasn.txt
