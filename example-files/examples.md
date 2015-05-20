__Example commands could be:__

```bash
/home/chrishah/Dropbox/Github/Docker/metab-box/scripts/metabbox-test/fix_gb.py -R REFfile -Q QUERYfile_full --taxids --blast --fasta --phyloplace --merge --product_length 400 --clust_cov 3 --min_ident 0.95 -E --trim_minlength 100 --n_threads 5 --gb_out ../references/CytB_cleaned_05_2015_checked.gb --rec_check

/home/chrishah/Dropbox/Github/Docker/metab-box/scripts/metabbox-test/fix_gb.py -Q QUERYfile_all -R REFfile --seqinfo --fasta --taxids --blast -m 12S -n 4 --min_ident 1 --clust_cov 3 --trim_adapter adapters_rc.fasta --trim_minlength 90 -E --merge --product_length 110 --rec_check
```


__for every query ID the pipeline will create a result directory__
example for sample `A01`

```bash
ls -hlrt A01/
total 127M
-rw-rw-r-- 1 chrishah chrishah 1.5M May  8 11:05 A01_forward.paired.fastq.gz
-rw-rw-r-- 1 chrishah chrishah 209K May  8 11:05 A01_forward.singletons.fastq.gz
-rw-rw-r-- 1 chrishah chrishah  40K May  8 11:05 A01_reverse.singletons.fastq.gz
-rw-rw-r-- 1 chrishah chrishah 1.5M May  8 11:05 A01_reverse.paired.fastq.gz
-rw-r--r-- 1 chrishah chrishah 760K May  8 11:05 A01.notCombined_2.fastq.gz
-rw-r--r-- 1 chrishah chrishah 771K May  8 11:05 A01.notCombined_1.fastq.gz
-rw-r--r-- 1 chrishah chrishah 927K May  8 11:05 A01.extendedFrags.fastq.gz
-rw-rw-r-- 1 chrishah chrishah 1.5K May  8 11:05 A01.histogram
-rw-rw-r-- 1 chrishah chrishah  184 May  8 11:05 A01.hist
-rw-rw-r-- 1 chrishah chrishah  16M May  8 11:05 A01_trimmed.fasta
-rw-rw-r-- 1 chrishah chrishah 6.7M May  8 11:05 A01.uc
-rw-rw-r-- 1 chrishah chrishah 896K May  8 11:05 A01_centroids_backup.fasta
-rw-rw-r-- 1 chrishah chrishah 169K May  8 11:05 A01_centroids.fasta
-rw-rw-r-- 1 chrishah chrishah  70M May  8 11:05 blastn.standard.out
-rw-rw-r-- 1 chrishah chrishah  29M May  8 11:05 marker_A01_blastn.out.xml
-rw-rw-r-- 1 chrishah chrishah 5.9K May  8 11:05 Gasterosteus_aculeatus.fasta
-rw-rw-r-- 1 chrishah chrishah  13K May  8 11:05 Abramis_brama.fasta
-rw-rw-r-- 1 chrishah chrishah  955 May  8 11:05 Salmo_salar.fasta
-rw-rw-r-- 1 chrishah chrishah 3.2K May  8 11:05 Rutilus_rutilus.fasta
-rw-rw-r-- 1 chrishah chrishah  84K May  8 11:05 Perca_fluviatilis.fasta
-rw-rw-r-- 1 chrishah chrishah 5.7K May  8 11:05 Salmo_trutta.fasta
-rw-rw-r-- 1 chrishah chrishah 2.0K May  8 11:05 Salmo.fasta
-rw-rw-r-- 1 chrishah chrishah 3.2K May  8 11:05 Percinae.fasta
-rw-rw-r-- 1 chrishah chrishah 1.1K May  8 11:05 Percidae.fasta
-rw-rw-r-- 1 chrishah chrishah 2.7K May  8 11:05 Cyprinidae.fasta
-rw-rw-r-- 1 chrishah chrishah  52K May  8 11:05 no_hit.fasta
-rw-rw-r-- 1 chrishah chrishah  457 May  8 11:05 A01-results.txt
```

__example results summary for sample A01__

```bash
cat A01-results.txt
```
```bash
The sample A01 contained 49510 valid query sequences:

species: 35343 (71.39 %)
        Abramis brama: 5865 (11.85 %)
        Gasterosteus aculeatus: 1202 (2.43 %)
        Perca fluviatilis: 26527 (53.58 %)
        Rutilus rutilus: 360 (0.73 %)
        Salmo salar: 19 (0.04 %)
        Salmo trutta: 1370 (2.77 %)
genus: 401 (0.81 %)
        Salmo: 401 (0.81 %)
subfamily: 770 (1.56 %)
        Percinae: 770 (1.56 %)
family: 1051 (2.12 %)
        Cyprinidae: 592 (1.20 %)
        Percidae: 459 (0.93 %)
nohit: 11528 (23.28 %)


Fri May  8 10:52:39 2015
```

__The overall log for a run could look as follows__

```bash
cat log
```
```bash
/home/chrishah/Dropbox/Github/Docker/metab-box/scripts/metabbox-test/fix_gb.py -R REFfile -Q QUERYfile_full --taxids --blast --fasta --phyloplace --merge --product_length 400 --clust_cov 3 --min_ident 0.95 -E --trim_minlength 100 --n_threads 5 --gb_out ../references/CytB_cleaned_05_2015_checked.gb --rec_check

processing /media/chrishah/STORAGE/eDNA/metabeatbox-analyses/CytB/analyses_05_2015/references/CytB_cleaned_05_2015.gb (containing 2155 records)

total number of records: 2155

write out taxids to taxids.txt

running taxit to generate reduced taxonomy table
taxit taxtable -d /home/chrishah/src/taxtastic/taxonomy_db/taxonomy.db -t taxids.txt
write out reference sequences to refs.fasta

building blast db for marker marker



Building a new DB, current time: 05/08/2015 11:05:02
New DB name:   marker_blast_db
New DB title:  refs.fasta
Sequence type: Nucleotide
Keep Linkouts: T
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 2155 sequences in 0.112985 seconds.


Fri May  8 11:05:02 2015


processing query ID: A01
###############


trimming PE reads with trimmomatic
java -jar /home/chrishah/src/Trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 5 -phred33 /home/chrishah/Dropbox/BaseSpace/Fish-Cyt-B-test-21170166/A01D-23323211/Data/Intensities/BaseCalls/A01D_S1_L001_R1_001.fastq.gz /home/chrishah/Dropbox/BaseSpace/Fish-Cyt-B-test-21170166/A01D-23323211/Data/Intensities/BaseCalls/A01D_S1_L001_R2_001.fastq.gz A01_forward.paired.fastq.gz A01_forward.singletons.fastq.gz A01_reverse.paired.fastq.gz A01_reverse.singletons.fastq.gz ILLUMINACLIP:None:5:5:5 TRAILING:30 LEADING:30 SLIDINGWINDOW:5:30 MINLEN:100

TrimmomaticPE: Started with arguments: -threads 5 -phred33 /home/chrishah/Dropbox/BaseSpace/Fish-Cyt-B-test-21170166/A01D-23323211/Data/Intensities/BaseCalls/A01D_S1_L001_R1_001.fastq.gz /home/chrishah/Dropbox/BaseSpace/Fish-Cyt-B-test-21170166/A01D-23323211/Data/Intensities/BaseCalls/A01D_S1_L001_R2_001.fastq.gz A01_forward.paired.fastq.gz A01_forward.singletons.fastq.gz A01_reverse.paired.fastq.gz A01_reverse.singletons.fastq.gz ILLUMINACLIP:None:5:5:5 TRAILING:30 LEADING:30 SLIDINGWINDOW:5:30 MINLEN:100
May 08, 2015 11:05:02 AM org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer makeIlluminaClippingTrimmer
SEVERE: null
java.io.FileNotFoundException: /media/chrishah/STORAGE/eDNA/metabeatbox-analyses/CytB/analyses_05_2015/trim_0.95_min_clust_3/A01/None (No such file or directory)
	at java.io.FileInputStream.open(Native Method)
	at java.io.FileInputStream.<init>(FileInputStream.java:146)
	at org.usadellab.trimmomatic.fasta.FastaParser.parse(FastaParser.java:54)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.loadSequences(IlluminaClippingTrimmer.java:107)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:70)
	at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:27)
	at org.usadellab.trimmomatic.TrimmomaticPE.run(TrimmomaticPE.java:495)
	at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:35)

Input Read Pairs: 39025 Both Surviving: 31720 (81.28%) Forward Only Surviving: 3386 (8.68%) Reverse Only Surviving: 701 (1.80%) Dropped: 3218 (8.25%)
TrimmomaticPE: Completed successfully

merging paired-end reads with flash

flash A01_forward.paired.fastq.gz A01_reverse.paired.fastq.gz -M 400 -t 5 -p 33 -o A01 -z
[FLASH] Starting FLASH v1.2.11
[FLASH] Fast Length Adjustment of SHort reads
[FLASH]  
[FLASH] Input files:
[FLASH]     A01_forward.paired.fastq.gz
[FLASH]     A01_reverse.paired.fastq.gz
[FLASH]  
[FLASH] Output files:
[FLASH]     ./A01.extendedFrags.fastq.gz
[FLASH]     ./A01.notCombined_1.fastq.gz
[FLASH]     ./A01.notCombined_2.fastq.gz
[FLASH]     ./A01.hist
[FLASH]     ./A01.histogram
[FLASH]  
[FLASH] Parameters:
[FLASH]     Min overlap:           10
[FLASH]     Max overlap:           400
[FLASH]     Max mismatch density:  0.250000
[FLASH]     Allow "outie" pairs:   false
[FLASH]     Cap mismatch quals:    false
[FLASH]     Combiner threads:      5
[FLASH]     Input format:          FASTQ, phred_offset=33
[FLASH]     Output format:         FASTQ, phred_offset=33, gzip
[FLASH]  
[FLASH] Starting reader and writer threads
[FLASH] Starting 5 combiner threads
[FLASH] Processed 25000 read pairs
[FLASH] Processed 31720 read pairs
[FLASH]  
[FLASH] Read combination statistics:
[FLASH]     Total pairs:      31720
[FLASH]     Combined pairs:   15873
[FLASH]     Uncombined pairs: 15847
[FLASH]     Percent combined: 50.04%
[FLASH]  
[FLASH] Writing histogram files.
[FLASH]  
[FLASH] FLASH v1.2.11 complete!
[FLASH] 0.286 seconds elapsed


zcat A01.notCombined_2.fastq.gz A01_reverse.singletons.fastq.gz | fastq_to_fasta | fastx_reverse_complement > temp2.fasta
zcat A01.extendedFrags.fastq.gz A01.notCombined_1.fastq.gz A01_forward.singletons.fastq.gz | fastq_to_fasta > temp1.fasta
cat temp1.fasta temp2.fasta > temp_trimmed.fasta

clustering using vsearch
vsearch --cluster_fast A01_trimmed.fasta --id 1.00 --threads 5 --centroids A01_centroids.fasta --uc A01.uc
vsearch v1.1.0_linux_x86_64, 31.4GB RAM, 12 cores
https://github.com/torognes/vsearch


vsearch identified 2691 clusters (clustering threshold 1.00) - 547 clusters (minimum of 3 records per cluster) are used in subsequent analyses

running blast search against local database ../marker_blast_db
blastn -query A01_centroids.fasta -db ../marker_blast_db -evalue 0.001 -outfmt 5 -out marker_A01_blastn.out.xml -num_threads 5 -max_target_seqs 50

query: M02538:27:000000000-ACTEL:1:1103:23826:7628_1:N:0:1
no hit - done

query: M02538:27:000000000-ACTEL:1:2110:14939:16749_1:N:0:1
no hit - done

query: M02538:27:000000000-ACTEL:1:1101:11701:5558_1:N:0:1
1 top90 match - done:
Perca fluviatilis
number of amgibuous hits: 0

query: M02538:27:000000000-ACTEL:1:1115:27474:20508_2:N:0:1
1 top90 match - done:
Salmo trutta
number of amgibuous hits: 0

query: M02538:27:000000000-ACTEL:1:1111:19695:17761_1:N:0:1
1 perfect match - done:
Salmo salar

query: M02538:27:000000000-ACTEL:1:1102:24675:7207_1:N:0:1
1 top90 match - done:
Gasterosteus aculeatus
number of amgibuous hits: 0

query: M02538:27:000000000-ACTEL:1:2103:17620:23801_2:N:0:1
1 top90 match - done:
Perca fluviatilis
number of amgibuous hits: 0

query: M02538:27:000000000-ACTEL:1:2119:28559:16950_2:N:0:1
number of amgibuous hits: 2
assigned to LCA - done
Cyprinidae

query: M02538:27:000000000-ACTEL:1:1104:28451:7814_1:N:0:1
number of amgibuous hits: 3
assigned to LCA - done
Percidae

query: M02538:27:000000000-ACTEL:1:1119:15857:21876_2:N:0:1
no hit - done

number of clusters processed: 547 (of 547)

The sample A01 contained 49510 valid query sequences:

species: 35343 (71.39 %)
        Abramis brama: 5865 (11.85 %)
        Gasterosteus aculeatus: 1202 (2.43 %)
        Perca fluviatilis: 26527 (53.58 %)
        Rutilus rutilus: 360 (0.73 %)
        Salmo salar: 19 (0.04 %)
        Salmo trutta: 1370 (2.77 %)
genus: 401 (0.81 %)
        Salmo: 401 (0.81 %)
subfamily: 770 (1.56 %)
        Percinae: 770 (1.56 %)
family: 1051 (2.12 %)
        Cyprinidae: 592 (1.20 %)
        Percidae: 459 (0.93 %)
nohit: 11528 (23.28 %)


Fri May  8 11:05:48 2015

```
