for i in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
    do python3 $HOME/lib/SIEVEKmers/KmerFeatures.py -f ../ubligase_examples_ids.fasta -R 9 -m gist -M reduced_alphabet_0 -k $i -o ubligase_k$i.red0;
done;
