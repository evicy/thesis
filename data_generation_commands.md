## Data generation commands
##### Sources:
[Reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2/)
[Available sample genome reads](https://www.ebi.ac.uk/ena/browser/view/PRJNA563564?fbclid=IwAR2X6qrgFk6290szGloBZXG2Bl_4hommUymZACM-ZSGq4OWwH2TZrjANxXw)

##### Tools used:
* [bwa](https://github.com/lh3/bwa)
* 
* [edso](https://github.com/webmasterar/edso)
---
### Get the vcf file

**Index the reference genome**
`bwa index genome_assemblies_genome_fasta.tar`

**Align and index every read (example with read [SRR10058833](https://www.ebi.ac.uk/ena/browser/view/SRR10058833))**
`bwa mem reference_genome.fna SRR10058833_1.fastq.gz SRR10058833_2.fastq.gz -R '@RG\tID:SRR10058833\tSM:SRR10058833\tPL:illumina\tLB:SRR10058833\tPU:unit' | samtools view -S -b - |  samtools sort -o SRR10058833.bam ; samtools index SRR10058833.bam`

**Create vcf file from the samples specified in bam_list.txt (bam_list.txt: line separated sample_i.bam)**
`freebayes -f  reference_genome.fna -L bam_list.txt > sample_i.vcf`

---
### Get the bubble graph
`./edso reference_genome.fna variants.vcf outfile.eds`

