# This file is not meant to be ran

# bwa can be used in the singularity wherever
# samtools can be used in the singularity wherever
# HLA-LA is stored in /usr/local/opt/hla-la/bin/HLA-LA

# directory permanent mount
/app/singularity <= program
/app/references <= references DB
/app/uploads <= user uploads

# bwa command
singularity exec -B /app/references:/tmp /app/singularity/bwa_0_7_18 bwa index /tmp/chr6.fa
singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa /tmp/uploads/16/e1.fq > /app/uploads/16/e1.sam
# samtools command
singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/uploads/16/e1.sam > /app/uploads/16/e1.bam
singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/uploads/16/e1.bam > /app/uploads/16/e1.bam.bai
# HLA-LA command
singularity exec -B /app:/tmp /app/singularity/hla-la_1_0_3 /usr/local/opt/hla-la/bin/HLA-LA --action testBinary
singularity exec -B /app:/tmp /app/singularity/hla-la_1_0_3 /usr/local/opt/hla-la/bin/HLA-LA --action prepareGraph --PRG_graph_dir /tmp/references/PRG_MHC_GRCh38_withIMGT

 singularity exec -B /app:/tmp /app/singularity/hla-la_1_0_3 /usr/local/opt/hla-la/src/HLA-LA.pl --BAM /tmp/uploads/16/e1.bam --graph /tmp/references/PRG_MHC_GRCh38_withIMGT --sampleID NA12878 --maxThreads 7 --workingDir /tmp