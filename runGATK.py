import subprocess
import argparse
import os

'''
Do prawidlowego dzialania skryptu wymagane jest aby:
- programy GATK i picard byly znajdowaly sie katalog wyzej
- w katalogu znajdowal sie katalog ref zawierajacy wymienione nizej pliki pochodzace z
  GATK bundle ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
        -Homo_sapiens_assembly38.fasta
        -Homo_sapiens_assembly38.dbsnp.vcf
        -hapmap_3.3.hg38.vcf
        -1000G_omni2.5.hg38.vcf
        -1000G_phase1.snps.high_confidence.hg38.vcf
        -Mills_and_1000G_gold_standard.indels.hg38.vcf
'''

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file",help="SAM file",type=str, required=True)
args = parser.parse_args()

fileName, fileExtension = os.path.splitext(args.file)
if fileExtension and fileExtension != '.sam':
    raise Exception('Error when parsing SAM file')

def runProcess(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

    if process.returncode == 1:
        exit(1)

#Add RG and coordinate sort
command = 'java -Xmx3g -jar ../picard.jar AddOrReplaceReadGroups I=' + fileName + '.sam O=' + fileName + '_mapped_rg.bam SORT_ORDER=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=I'
print command
runProcess(command)

#Mark duplicates
command = 'java -Xmx3g -jar ../picard.jar MarkDuplicates I=' + fileName + '_mapped_rg.bam O=' + fileName + '_marked_duplicates_rg.bam M=' + fileName + '_marked_dup_metrics_rg.txt CREATE_INDEX=true'
print command
runProcess(command)

#Create target list of intervals to be realigned
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref/Homo_sapiens_assembly38.fasta -I ' + fileName + '_marked_duplicates_rg.bam -o ' + fileName + '_forIndelRealigner.intervals.list'
print command
runProcess(command)

command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R ref/Homo_sapiens_assembly38.fasta -I ' + fileName + '_marked_duplicates_rg.bam -targetIntervals ' + fileName + '_forIndelRealigner.intervals.list -o ' + fileName + '_realignedBam.bam'
print command
runProcess(command)

#Base recalibration
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R ref/Homo_sapiens_assembly38.fasta -I ' + fileName + '_realignedBam.bam -knownSites ref/Homo_sapiens_assembly38.dbsnp.vcf -o ' + fileName + '_recal_data.table'
print command
runProcess(command)

#Apply the recalibration to sequence data
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads -R ref/Homo_sapiens_assembly38.fasta -I ' + fileName + '_realignedBam.bam -BQSR ' + fileName + '_recal_data.table -o ' + fileName + '_recal_reads.bam'
print command
runProcess(command)

#Variant calling
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -I ' + fileName + '_recal_reads.bam -R ref/Homo_sapiens_assembly38.fasta --output_mode EMIT_VARIANTS_ONLY -ploidy 2 -o ' + fileName + '_raw_variants.vcf'
print command
runProcess(command)

#dbSNP annotation
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantAnnotator -R ref/Homo_sapiens_assembly38.fasta --variant ' + fileName + '_raw_variants.vcf --dbsnp ref/Homo_sapiens_assembly38.dbsnp.vcf -L ' + fileName + '_raw_variants.vcf -o ' + fileName + '_annotated_raw_variants.vcf'
print command
runProcess(command)

#Variant recalibrator
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantRecalibrator -R ref/Homo_sapiens_assembly38.fasta -input ' + fileName + '_raw_variants.vcf  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/hapmap_3.3.hg38.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 ref/1000G_omni2.5.hg38.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/1000G_phase1.snps.high_confidence.hg38.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/Homo_sapiens_assembly38.dbsnp.vcf -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile ' + fileName + '_recalibrate_SNP.recal -tranchesFile ' + fileName + '_recalibrate_SNP.tranches'
print command
runProcess(command)

#Apply SNP Recalibration
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T ApplyRecalibration  -R ref/Homo_sapiens_assembly38.fasta -input ' + fileName + '_raw_variants.vcf  -mode SNP --ts_filter_level 99.0 -recalFile ' + fileName + '_recalibrate_SNP.recal -tranchesFile ' + fileName + '_recalibrate_SNP.tranches -o ' + fileName + '_recalibrated_snps_raw_indels.vcf'
print command
runProcess(command)

#Build the Indel recalibration model
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantRecalibrator  -R ref/Homo_sapiens_assembly38.fasta -input ' + fileName + '_recalibrated_snps_raw_indels.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 ref/Mills_and_1000G_gold_standard.indels.hg38.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/Homo_sapiens_assembly38.dbsnp.vcf -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile ' + fileName + '_recalibrate_INDEL.recal -tranchesFile ' + fileName + '_recalibrate_INDEL.tranches'
print command
runProcess(command)

#Apply the desired level of recalibration to the Indels in the call set
command = 'java -Xmx3g -jar ../GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T ApplyRecalibration  -R ref/Homo_sapiens_assembly38.fasta -input ' + fileName + '_recalibrated_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile ' + fileName + '_recalibrate_INDEL.recal -tranchesFile ' + fileName + '_recalibrate_INDEL.tranches -o ' + fileName + '_recalibrated_variants.vcf'
print command
runProcess(command)

print 'ok'
