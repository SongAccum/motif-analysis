# code used for motif analysis
## motif.sh
* get `data needed` for subsequent analysis(sequence & structure motif) pipeline
* only need raw differentially expressed gene lists
* change homefolder
## sequence_motif.sh
* perform `sequence motif` analysis pipeline
* split from motif.sh to run more efficiency
* meme de novo motif discovery
* ame motif enrichment
## structure_motif.sh
* perform `structure motif` analysis pipeline
* split from motif.sh to run `more efficiency`
* **RNAfold**
* **beam**
* **weblogo**
* **varna**
## concat.r
* `concatenate sequences` of the same transcripts
* passing parameter `\*.fa`
* change the `source` file 
## bg_down.r
* get downstream `1000bp` as bg sequence
* passing parameters input:`\*.bed` output:`\*_bg.bed`
## visual_2D.py
* perform `VARNA: Visualization Applet for RNA` pipeine
