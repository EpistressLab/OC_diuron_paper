import os
import sys
import errno
import json

# def getFastqFileDict():
#     fileDict={}
#     for s in config["samples"]:
#         fileDict[s['sample']]=[s['files']['forward'],s['files']['reverse']]
#     return fileDict
def getFastqFileFW(wildcards):
    for s in config["samples"]:
        if (s["sample"]== wildcards.sampleName):
            return s["files"]["forward"]

def getFastqFileRV(wildcards):
    for s in config["samples"]:
        if (s["sample"]== wildcards.sampleName):
            return s["files"]["reverse"]

def getConditions():
	sampleList = list()
	for cond in config["replicates"].keys():
		for s in config["replicates"][cond]["samples"]:
			sd = dict()
			sd["samplename"] = s
			sd["filename"] = s
			sd["condition"] = cond
			sampleList.append(sd)
	return sampleList
####################	 CMD to launch pipeline   ############################
## snakemake --config projf="project.json" reff="references.json" conff="config.json" -nrp

########## loading configuration variable ##################
# configfile: config["reff"]
configfile: config["projf"]
configfile: config["conff"]


# forwardFiles = list()
# reverseFiles = list()
sample_all = list()
# fileDict={}
for s in config["samples"]:
    # s=config["samples"][i]
    sample_all.append(s["sample"])
    # fileDict.update( {s['sample']:[s['files']['forward'],s['files']['reverse']]} )
# fileDict=getFastqFileDict()

# sampleNames=["HOS_NI_2","HOS_shLacZ_4","HOS_shLacZ_5","HOS_shHSF1_482_7","HOS_shHSF1_482_8","HOS_shHSF1_483_11",]
OUTPUTDIR="Results_ncbi37.2"
GTF=config["reference"]["gtf"]

HS2_GENOME_INDEX = config["reference"]["hisat2_genome_index_base"]
FASTA = config["reference"]["fasta"]
SCRIPTPATH = "/home/rob/Documents/tools/RNASeq/RNAseq_quantif_pipeline/scripts"

BIOMART = config["reference"]["biomart"]
REFNAME = config["reference"]["name"]

LOGFCTHRESHOLD=1
FDRTHRESHOLD=0.05

rule all:
    input: 
        expand(OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_{read}.fastq.gz",sampleName=sample_all,read=["R1","R2"]),
        expand(OUTPUTDIR+"/DESEQ2/featureCounts/{sampleName}",sampleName=sample_all),
        expand(OUTPUTDIR+"/DESEQ2/counts_featC/{sampleName}",sampleName=sample_all),
        # #OUTPUTDIR+"/DESEQ2/results_FDRThresh5/NormalizedCountMatrix.txt",
        OUTPUTDIR+"/DESEQ2/DESEQ2_CONDITIONS_all.tab"
# rule filterPrinseq:
#     input: R1=getFastqFileFW,
#         R2=getFastqFileRV
#     output: R1=OUTPUTDIR+"/Samples/{sampleName}/PRINSEQ/{sampleName}_good_1.fastq",
#         R2=OUTPUTDIR+"/Samples/{sampleName}/PRINSEQ/{sampleName}_good_2.fastq"
#     params:	preGood = OUTPUTDIR+"/Samples/{sampleName}/PRINSEQ/{sampleName}_good",
#         preBad = OUTPUTDIR+"/Samples/{sampleName}/PRINSEQ/{sampleName}_bad",
#         log = OUTPUTDIR+"/Samples/{sampleName}/PRINSEQ/{sampleName}.log"
#     shell:"""
#         prinseq-lite.pl -fastq {input.R1} -fastq2 {input.R2} -out_good {params.preGood} -out_bad null -log {params.log} -no_qual_header -min_len 30 -max_len 60 -min_qual_mean 30 \
#         -ns_max_n 25 -graph_data {OUTPUTDIR}+"/Samples/{wildcards.sampleName}/PRINSEQ/{wildcards.sampleName}.gd -qual_no_scale -stats_all \
#         -derep 14 -derep_min 6 > OUTPUTDIR+"/Samples/{wildcards.sampleName}/PRINSEQ/{wildcards.sampleName}.txt
    # """
# rule filterFastxQual:
#     input: R1=getFastqFileFW,
#         R2=getFastqFileRV
#     output: R1=OUTPUTDIR+"/Samples/{sampleName}/FASTX/{sampleName}_R1.fastq.gz",
#         R2=OUTPUTDIR+"/Samples/{sampleName}/FASTX/{sampleName}_R2.fastq.gz"
#     shell:"""
#         fastq_quality_filter -q 30 -p 50 -z -i {input.R1} -o {output.R1} > {wildcards.sampleName}_R1.fastx_report &&
#         fastq_quality_filter -q 30 -p 50 -z -i {input.R2} -o {output.R2} > {wildcards.sampleName}_R2.fastx_report
#     """
rule filterFastp:
    input: R1=getFastqFileFW,
        R2=getFastqFileRV
    output: R1=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_R1.fastq.gz",
        R2=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_R2.fastq.gz",
        report=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_report.txt",
        html=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}.html"
    threads: 6
    shell:"""
        fastp -i {input.R1} -o {output.R1} -I {input.R2} -O {output.R2} --detect_adapter_for_pe -q 30 --html {output.html} \
        --thread {threads} > {output.report}
    """

rule align_PE_hisat2:
    input:
        R1=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_R1.fastq.gz",
        R2=OUTPUTDIR+"/Samples/{sampleName}/FASTP/{sampleName}_R2.fastq.gz"
    output: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.sam"
	threads: 12
	# conda:	"CONDA/rnaSeqQuantif.yml"
	shell:"""
        hisat2 --threads {threads} --rna-strandness RF --downstream-transcriptome-assembly --rg "SM:{wildcards.sampleName}" \
		--rg-id {wildcards.sampleName} --dta -x {HS2_GENOME_INDEX} -1 {input.R1} -2 {input.R2} \
        -S {OUTPUTDIR}/Samples/{wildcards.sampleName}/hisat2/{wildcards.sampleName}.sam --met-file \
        {OUTPUTDIR}/Samples/{wildcards.sampleName}/hisat2/{wildcards.sampleName}_metrics.txt \
        --summary-file {OUTPUTDIR}/Samples/{wildcards.sampleName}/hisat2/{wildcards.sampleName}_summary.txt
    """
rule sortSam:
	input: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.sam"
	output: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.poSrt.sam"
	threads: 6
	shell:"""
		samtools sort --threads {threads} -o {output} {input} 
		"""
# rule htseqCount:
#     input: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.poSrt.sam"
#     output: sam=OUTPUTDIR+"/Samples/{sampleName}/htseqCount2/{sampleName}.sam",
#         counts=OUTPUTDIR+"/DESEQ2/htseqCount2/{sampleName}"
#     shell:"""
#         htseq-count -o {output.sam} -r pos -s reverse {input} {GTF} > {output.counts}
#     """
# rule featureCounts:
#     input: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.poSrt.sam"
#     output: OUTPUTDIR+"/DESEQ2/featureCounts/{sampleName}",
#     threads: 6
#     shell:"""
#         featureCounts -a {GTF} -o {output} -O --nonOverlap 5 -s 2 -p -B -T {threads} --extraAttributes gene_name {input}
#     """

rule featureCounts2:
    input: OUTPUTDIR+"/Samples/{sampleName}/hisat2/{sampleName}.poSrt.sam"
    output: OUTPUTDIR+"/DESEQ2/featureCounts/{sampleName}",
    threads: 6
    shell:"""
        featureCounts -a {GTF} -o {output} -J -G {FASTA} -s 2 -p -B -T {threads} --extraAttributes gene_name {input}
    """

# rule prepareCountData:
#     input: OUTPUTDIR+"/DESEQ2/featureCounts/{sampleName}"
#     output: OUTPUTDIR+"/DESEQ2/counts_featC/{sampleName}"
#     shell:"""
#         awk '/ENSG/ {{print $1, $NF}}' {input} > {output}
#     """

rule prepareCountData2:
    input: OUTPUTDIR+"/DESEQ2/featureCounts/{sampleName}"
    output: OUTPUTDIR+"/DESEQ2/counts_featC/{sampleName}"
    shell:"""
        tail -n +3 {input} | awk '{{print $1, $NF}}'> {output}
    """
# rule prepareCountData:
#     input: OUTPUTDIR+"/DESEQ2/htseqCount2/{sampleName}"
#     output: OUTPUTDIR+"/DESEQ2/counts2/{sampleName}"
#     shell:"""
#         awk '{{/ENSG/ printf $1, $NF}}' {input} > {output}
#     """


# rule deseq2_conditions:
# 	output: tab = OUTPUTDIR+"/DESEQ2/DESEQ2_CONDITIONS.tab"
# 	run:
# 		condArray = getConditions()
# 		with open(output.tab, 'w') as condfile:
# 			condfile.write("samplename"+"\t"+"filename"+"\t"+"condition"+"\n")
# 			for s in condArray:
# 				condfile.write(s["samplename"]+"\t"+s["filename"]+"\t"+s["condition"]+"\n")

rule deseq2_conditions_2:
	output: tab = OUTPUTDIR+"/DESEQ2/DESEQ2_CONDITIONS_all.tab"
	run:
		condArray = getConditions()
		with open(output.tab, 'w') as condfile:
			condfile.write("samplename"+"\t"+"filename"+"\t"+"condition"+"\n")
			for s in condArray:
				condfile.write(s["samplename"]+"\t"+s["filename"]+"\t"+s["condition"]+"\n")

# rule deseq2:
# 	input:
# 		conditions=OUTPUTDIR+"/DESEQ2/DESEQ2_CONDITIONS.tab",
# 		counts=expand(OUTPUTDIR+"/DESEQ2/counts_featC/{sampleName}",sampleName=sample_all),
# 		corrAnnotations=(SCRIPTPATH+"/corresIDorg.txt")
# 	output: expand(OUTPUTDIR+"/DESEQ2/results_FDRThresh5/{suffixe}",suffixe=["NormalizedCountMatrix.txt","NormalizedCountMatrixFiltered.txt","PCAplot.png","sampletosampledistance.jpeg"]),
# #	conda: 	"CONDA/deseq2.yml"
# 	shell:	"""
# 		cat {SCRIPTPATH}/run_deseq2.R | R --slave --args {input.conditions} {OUTPUTDIR}/DESEQ2/counts_featC {OUTPUTDIR}/DESEQ2/results_FDRThresh5 {BIOMART} {REFNAME} {input.corrAnnotations} {LOGFCTHRESHOLD} {FDRTHRESHOLD}
# 		"""

# rule deseq2_2:
# 	input:
# 		conditions=OUTPUTDIR+"/DESEQ2/DESEQ2_CONDITIONS_all.tab",
# 		counts=expand(OUTPUTDIR+"/DESEQ2/counts_featC2/{sampleName}",sampleName=sample_all),
# 		corrAnnotations=(SCRIPTPATH+"/corresIDorg.txt")
# 	output: expand(OUTPUTDIR+"/DESEQ2/results_all_FDRThresh5/{suffixe}",suffixe=["NormalizedCountMatrix.txt","NormalizedCountMatrixFiltered.txt","PCAplot.png","sampletosampledistance.jpeg"]),
# #	conda: 	"CONDA/deseq2.yml"
# 	shell:	"""
# 		cat {SCRIPTPATH}/run_deseq2.R | R --slave --args {input.conditions} {OUTPUTDIR}/DESEQ2/counts_featC2 {OUTPUTDIR}/DESEQ2/results_all_FDRThresh5 {BIOMART} {REFNAME} {input.corrAnnotations} {LOGFCTHRESHOLD} {FDRTHRESHOLD}
# 		"""