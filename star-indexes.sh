#PBS -l select=1:mem=60gb:ncpus=6
#PBS -l walltime=05:00:00
#PBS -N starIndex

module load star
DIR=/rds/general/user/hm1412/projects/cebolalab_liver_regulomes/ephemeral/reference-genomes/GRCh38

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $DIR --genomeFastaFiles $DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

mv * $DIR
