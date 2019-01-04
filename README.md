# motif analysis
This is a document for motif related analysis.
## get 3'-UTR and promoter sequence
### install R package GenomicFeatures and biozhuoer tools
* **GenomicFeatures** package used to extract needed sequence
* **biozhuoer** tools used to concat sequences of the same 3’ UTR or promoter
```r
source("http://www.bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
library(GenomicFeatures)

if (!('devtools' %in% .packages(T))) install.packages('devtools');
devtools::install_github('dongzhuoer/biozhuoer');
```
### generate txdb object
There are many functions for us to get genme annotation file:
* makeTxDbFromBiomart  
* makeTxDbFromEnsembl  
* makeTxDbFromGFF      
* makeTxDbFromGRanges  
* makeTxDbFromUCSC
Here we use makeTxDbFromGFF for the existing of gtf file.
```r
gtf_file="/BioII/lulab_b/songyabing/genome/gencode.v27.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
```
### get 3'UTR & 5'UTR site range
```r
utr5p = fiveUTRsByTranscript(txdb, use.names=T)
utr3p = threeUTRsByTranscript(txdb, use.names=T)

utr3p.df=as.data.frame(utr3p)
utr5p.df=as.data.frame(utr5p)

write.table(utr3p.df, "utr3p.info", row.names=FALSE, sep='\t',quote=FALSE )
write.table(utr5p.df, "utr5p.info", row.names=FALSE, sep='\t' ,quote=FALSE)
```
### get promoter site range
```r
promoter=promoters(txdb)
promoter.df=as.data.frame(promoter)
write.table(promoter.df, "promoter.info", row.names=FALSE, sep='\t' ,quote=FALSE)
```
## intersect with diff_expr_gene
intersect with diff_exp_genes & get interested genes' info      
### interested 3'UTR
* 1st column:chr
* 2nd column:gene
* 3rd column:transprict
* 4th column:exon_name 
* 5th column:start
* 6th column:end
* 7th column:width(seq length)
* 8th column:strand
* 9th column:exon_id
* 10th column:exon_rank
```bash
sort -t $'\t' -k 2 utr3p.info|join -o 1.3 2.1 1.2 1.9 1.4 1.5 1.6 1.7 1.8 1.10 -t $'\t' -1 2 -2 2 - \
  <(cut -f 1 SC2_SF2.ct.dn.1_0.01.protein_coding |sort |join -t $'\t' -1 1 -2 1 - <(sort -t $'\t' -k 1 \
  <(grep -o -P -e "gene_id.*; transcript_id.*?;" gencode.v27.annotation.gtf |sort \
  |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//'))|sort -t $'\t' -k 2  ) |\
  sort -t $'\t' -k 1 >interested_three_prime_UTR.info
```
### interested promoter
* 1st column:chr
* 2nd column:gene
* 3rd column:transprict_name
* 4th column:start
* 5th column:end
* 6th column:width(seq length)
* 7th column:strand
* 8th column:transprict_id
```bash
sort -t $'\t' -k 7 promoter.info|join -o 1.1 2.1 1.7 1.2  1.3 1.4 1.5 1.6 -t $'\t' -1 7 -2 2 - \
  <(cut -f 1 SC2_SF2.ct.dn.1_0.01.protein_coding |sort |join -t $'\t' -1 1 -2 1 - \
  <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" gencode.v27.annotation.gtf |sort \
  |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//' ))|sort -t $'\t' -k 2  ) |\
  sort -t $'\t' -k 1 >interested_promoter.info
```
## convert to BED format
### 3'UTR bed info 
* 1st column:chr
* 2nd column:start
* 3rd column:end
* 4th column:gene
* 5th column:transprict
* 6th column:strand

**change for subsequent concat seq where 4th column: transcript 5th column: gene**
```bash
cat interested_three_prime_UTR.info | \
  awk '{print $1 "\t" $5-1 "\t" $6 "\t" $3 "\t" $2 "\t" $8}' | \
  sort -u  > interested_three_prime_UTR.bed
```
### promoter bed info 
* 1st column:chr 
* 2th column:start
* 3th column:end
* 4nd column:gene
* 5rd column:transprict_name
* 6th column:strand

**change for subsequent concat seq where 4th column: transcript 5th column: gene**
```bash
cat interested_promoter.info | \
  awk '{print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" $2 "\t" $7}' | \
  sort -u  > interested_promoter.bed
```
## get genome sequence
### get 3'UTR related genome sequence
* -s: Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. 
* -name:	Use the “name” column in the BED file for the FASTA headers in the output FASTA file.
* -fi: input FASTA
* -fo: Specify an output file name. By default, output goes to stdout.
```{bash,eval=FALSE}
bedtools getfasta -s -name -fi GRCh38.p10.genome.fa \
  -bed interested_three_prime_UTR.bed -fo interested_three_prime_UTR.fa
```
### concatenate sequences of the same 3’ UTR
```r
concatenate_seq <- function(fasta_file) {
    biozhuoer::read_fasta(fasta_file) %>% 
        dplyr::mutate(name = stringr::str_extract(name, 'ENST[\\d\\.]+')) %>% 
        dplyr::group_by(name) %>% dplyr::summarise(seq = paste0(seq, collapse = '')) %>% 
        biozhuoer::write_fasta(fasta_file)
}
concatenate_seq('interested_three_prime_UTR.fa')
```

### get promoter related genome sequence
```bash
bedtools getfasta -s -name -fi GRCh38.p10.genome.fa \
  -bed interested_promoter.bed -fo interested_promoter.fa
```

### concatenate sequences of the same promoter
For the purpose of concat transcript seq with same promoter ,we should change the columns order of the bed file **(with the 4rd column is transcript name)**.
```{r,eval=FALSE}
concatenate_seq <- function(fasta_file) {
    biozhuoer::read_fasta(fasta_file) %>% 
        dplyr::mutate(name = stringr::str_extract(name, 'ENST[\\d\\.]+')) %>% 
        dplyr::group_by(name) %>% dplyr::summarise(seq = paste0(seq, collapse = '')) %>% 
        biozhuoer::write_fasta(fasta_file)
}
concatenate_seq('interested_promoter.fa')
```

## motif enrichment with ame in meme suite

### generate random sequence
there are three mothods to get random sequence: 

* shuffle the input sequence 
* downsteam 1000bp  
* bedtools shuffle

#### shuffle the input sequence
```bash
fasta-shuffle-letters interested_three_prime_UTR.fa interested_three_prime_UTR.control

fasta-shuffle-letters interested_promoter.fa interested_promoter.control
```

#### downstream 1000bp as bg 
* https://dongzhuoer.github.io/diff_exp_2018_zhuoer/motif.html 
```r
slide <- function(input_bed, output_bed, n = 1000) {
    col_names <- c('chr', 'start', 'end', 'name', 'score', 'strand');

    original <- readr::read_tsv(input_bed, col_names) %>% 
        dplyr::group_by_at(-2:-3) %>% 
        dplyr::summarise(length = sum(end - start), end = max(end)) %>% 
        dplyr::ungroup()
    
    if (n > 0) {
       slide <- original %>% dplyr::mutate(start = end + n, end = start + length)
    } else {
       slide <- original %>% dplyr::mutate(end = start + n, start = end - length)
    }
    
    slide %>% dplyr::select(chr, start, end, name, score, strand) %>% 
        readr::write_tsv(output_bed, col_names = F)
}
slide('interested_three_prime_UTR.bed', 'interested_three_prime_UTR_downstream.bed')
slide('interested_promoter.bed', 'interested_promoter_downstream.bed')
```
**repeat  get promoter and get 3'UTR section**

#### bedtools shuffle
```bash
bedtools shuffle -i interested_three_prime_UTR.bed \
-g GRCh38.p10.genome.size >interested_three_prime_UTR_btools.bed

bedtools shuffle -i interested_promoter.bed \
-g GRCh38.p10.genome.size >interested_promoter_btools.bed
```
**repeat  get promoter and get 3'UTR section**

### get de novo motif
* -nmotifs: the number of motifs to be reported, default value is one, remember to change
* -evt: MEME will stop searching for motifs if the last motif found has an E-value > ev.
* -minw： Search for motifs with a width ≥ minw. default:8
* -minw: Search for motifs with a width ≤ maxw. default:50
* -dna: DNA sequence 	alphabet
* -mod: used to describe the distribution of motif sites. 
```bash
meme -dna -minw 4 -minw 12 -oc UTR_de_novo -mod zoops -nmotifs 5 interested_three_prime_UTR.fa

#run with a bug:Dataset too large (> 100000).  Rerun with larger -maxsize.
#version too old
#meme -dna -oc /BioII/lulab_b/songyabing/motif_analysis/motif_db/promoter_de_novo -mod zoops -nmotifs 30 /BioII/lulab_b/songyabing/motif_analysis/data/interested_promoter.fa

meme -dna -maxsize 1000000 \
-minw 4 -minw 12 \
-oc promoter_de_novo \
-mod zoops -nmotifs 5 \
interested_promoter.fa
```
### motif enrichment with ame
* **plus de novo motif file by meme**
* **do not need to convert**

* plus de novo motif file by meme
```bash
ame --control interested_three_prime_UTR.control  \
--oc UTR_output interested_three_prime_UTR.fa \
Homo_sapiens.meme Ray2013_rbp_Homo_sapiens.meme UTR_de_novo/meme.txt
```
* plus de novo motif file by meme
```bash
ame --control interested_promoter.control \
--oc promoter_output interested_promoter.fa \
JASPAR2018_CORE_vertebrates_non-redundant.meme Homo_sapiens.meme HOCOMOCOv11_core_HUMAN_mono_meme_format.meme promoter_de_novo/meme.txt
```
