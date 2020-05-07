import os
from fnmatch import fnmatch
import pandas as pd
from subprocess import call, check_output
from Bio import SeqIO

# Inputs and required files for pipeline
in_dir=config["in_dir"]
out_dir=config["out_dir"]
ref_sequence=config["ref_sequence"]
itype=config["barcode_check"]["type"]
limit=config["barcode_check"]["limit"]
barcodes=config["barcode_check"]["barcodes"]
barcode_regions=config["barcode_check"]["barcode_regions"]

# Find all gzipped files in input directory and add to sample dictionary.
# Within dictionary, group files based on their sample. I.e. group first and second reads.
SAMPLES = dict()
for file in os.listdir( in_dir ):
    if fnmatch( file, "*.gz" ):
        sample_name = file.split( "_" )[0]
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = [file]
        else:
            SAMPLES[sample_name].append( file )

# Load sample dataset, and generate final file names according to metadata.
# Will fit format: W###_Collection-date_Country_State_County_Latitude_Longitude
california = pd.read_csv( os.path.join( config["res"], "california_samples.csv" ), usecols=["Scripps_ID", "Collection date", "Country", "State", "County", "Latitude", "Longitude"] )
california["Collection date"] = pd.to_datetime( california["Collection date"] )
california["Collection date"] = california["Collection date"].dt.strftime( "%Y-%m-%d" )
california["Latitude"] = pd.to_numeric( california["Latitude"], errors="coerce" )
california["Latitude"] = california["Latitude"].round( 3 ).apply( lambda x: "{:.3f}".format( x ) )
california["Longitude"] = pd.to_numeric( california["Longitude"], errors="coerce" )
california["Longitude"] = california["Longitude"].round( 3 ).apply( lambda x: "{:.3f}".format( x ) )
california["File"] = california.apply( lambda x: ("_".join( [str( i ) for i in x] ) ).strip(), axis=1, raw=True )
california.index = california["Scripps_ID"]

# Assign final file names to each sample.
final_file = dict()
for i in SAMPLES.keys():
    try:
        final_file[i] = california.loc[i,"File"]
    except KeyError:
        final_file[i] = i

rule all:
    """ Output rule; specifies files we want to generate

        Outputs:
        - {sample}.fa                           : Consensus sequence for each sample
        - {sample}_bs.bam                       : Alignment of unmapped reads to reference barcode sequences
        - {full_name}.fasta                     : Consensus sequence for each sample, renamed to usable format
        - {sample}.coverage.png                 : Plot showing read coverage at each base of reference genome
        - {sample}.aligned.sorted.bam           : All reads aligned to reference genome
        - {sample}.trimmed.aligned.sorted.bam   : Reads aligned to reference genome and quality trimmed
        - alignment_statistics.csv              : Records a number of statistics for each sample
        - barcode_statistics.csv                : Records the number of reads which align to each barcode in barcode reference"""

    input:
        expand( "{out_dir}/_consensus/{sample}.fa", out_dir = out_dir, sample = SAMPLES ),
        expand( "{out_dir}/_barcode/{sample}_bs.bam", out_dir=out_dir, sample=SAMPLES ),
        expand( "{out_dir}/_final/{full_name}.fasta", out_dir=out_dir, full_name=final_file.values() ),
        expand( "{out_dir}/_coverage/{sample}.coverage.png", out_dir = out_dir, sample = SAMPLES ),
        expand( "{out_dir}/_untrimmed/{sample}.aligned.sorted.bam", out_dir = out_dir, sample = SAMPLES ),
        expand( "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam", out_dir = out_dir, sample = SAMPLES ),
        os.path.join( out_dir, "alignment_statistics.csv" ),
        os.path.join( out_dir, "barcode_statistics.csv" )
    
rule generate_statistics:
    """ Calculates a number of statistics for each sample. Calculates alignment statistics like number of reads, average per base coverage, and reference genome coverage. Also colates barcode statistics into a single data set.

        Inputs:
        - {sample}.trimmed.aligned.sorted.bam   : Reads aligned to reference genome and quality trimmed

        Outputs:
        - alignment_statistics.csv  : Records a number of statistics for each sample
        - barcode_satistics         : Records the number of reads which align to each barcode in barcode reference"""
    input:
        align = expand( "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam", out_dir=out_dir, sample=SAMPLES ),
        barcodes = expand( "{out_dir}/_barcode/{sample}.csv", out_dir=out_dir, sample=SAMPLES )
    output:
        alignment_stats = "{out_dir}/alignment_statistics.csv",
        barcode_stats = "{out_dir}/barcode_statistics.csv"
    run:
        # Alignment Statistics
        data_dict = { "Sample" : list(), "Total_Reads" : list(), "Percent_Coverage" : list(), "Mean_Coverage" : list() }
        for i in input.align:
            sample = os.path.splitext( os.path.basename( i ) )[0]
            data_dict["Sample"].append( sample )
            data_dict["Total_Reads"].append( check_output( ["samtools", "view", "-c", i] ).strip() )
            
            coverage = check_output( ["samtools", "depth", "-aa", "-d", "0", i] )

            bps = 0
            mean_coverage = 0
            percent_coverage = 0
            low_coverage = 0
            count = 0

            for line in coverage.decode().split( "\n" ):
                if line:
                    line_split = line.split( "\t" )
                    if 96 < int( line_split[1] ) < 10399:
                        count += 1
                        cov = int( line_split[2] )
                        mean_coverage += cov

                        if cov < 10:
                            low_coverage += 1
            
            data_dict["Mean_Coverage"].append( mean_coverage / count )
            data_dict["Percent_Coverage"].append( ( float( count - low_coverage ) / count ) * 100 )
            
        data_df = pd.DataFrame( data_dict )
        data_df.to_csv( str( output.alignment_stats ) )

        # Barcode Stats
        bc_list = list()
        for i in input.barcodes:
            temp_bc = pd.read_csv( i )
            temp_bc = temp_bc.set_index( "Barcode" )
            bc_list.append( temp_bc )
        bc = pd.concat( bc_list, axis=1 )
        bc.to_csv( output.barcode_stats )
        
rule cleanup_consensus:
    """ Renames consensus sequences so that they match a useful format. The default output of iVar isn't the best if multiple samples are run in parallel, so I rename them here. 

        Input:
        - {sample}.fa   : Consensus sequence directly from iVar consensus

        Output:
        - {full_name}.fasta : Renamed consensus according to format W###_Collection-date_Country_State_County_Latitude_Longitude """

    input:
        expand( "{out_dir}/_consensus/{sample}.fa", out_dir=out_dir, sample=SAMPLES )
    output:
        expand( "{out_dir}/_final/{full_name}.fasta", out_dir=out_dir, full_name=final_file.values() )
    run: 
        for i in input:
            sample = os.path.splitext( os.path.basename( i ) )[0]
            record = SeqIO.read( i, "fasta" )
            record.id = final_file[sample]
            record.description = ""
            SeqIO.write( record, os.path.join( out_dir, "_final/{}.fasta".format( record.id ) ), "fasta" )

rule generate_consensus:
    """ Generates a consensus sequence from an alignment.

        Input:
        - {sample}.trimmed.aligned.sorted.bam   : Reads aligned to reference sequence and trimmed for quality and primer sequences

        Output:
        - {sample}.fa   : Consensus sequence in fasta format"""

    input:
        "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam"
    output:
        "{out_dir}/_consensus/{sample}.fa"
    shell:
        "samtools mpileup -A -d 300000 -Q 0 -F 0 {input} | ivar consensus -m 10 -n N -p {output}"

rule generate_coverage_plot:
    input:
        "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam"
    output:
        "{out_dir}/_coverage/{sample}.coverage.csv",
        "{out_dir}/_coverage/{sample}.coverage.png"
    shell:
        "samtools depth -aa -d 0 {input} > {output[0]} &&"
        "python {config[scripts]}/graph_coverage.py {output[0]} {output[1]}"

rule trim_primer_quality:
    input: 
        "{out_dir}/_untrimmed/{sample}.aligned.sorted.bam" 
    output:
        "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam"
    shell:
        "ivar trim -b {config[bed_file]} -p {out_dir}/_aligned_bams/{wildcards.sample}.trimmed.aligned.sorted -i {input} &&"
        "samtools sort -o {output} {output} &&"
        "samtools index {output}"

rule barcode_check:
    input:
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][0],
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][1]
    output:
        unmapped = temp( "{out_dir}/_barcode/{sample}.bam" ),
        unmapped_r1 = temp( "{out_dir}/_barcode/{sample}_r1.fastq" ),
        unmapped_r2 = temp( "{out_dir}/_barcode/{sample}_r2.fastq" ),
        barcode_bam = "{out_dir}/_barcode/{sample}_bs.bam",
        barcode_stats = "{out_dir}/_barcode/{sample}.csv"
    run:
        # Load required modules
        command = "module load bedtools && "
        if type == "amplicon":
            command += "module load pear && "

        # Align reads to west nile genome and extract unmapped reads
        command += "bwa mem {} {} {} | samtools view -b -f 12 -F 256 | samtools sort -o {} && ".format( ref_sequence, input[0], input[1], output.unmapped )
        command += "samtools index {} && ".format( output.unmapped )

        # Seperate read pairs into seperate files from alignment
        command += "bedtools bamtofastq -i {} -fq {} -fq2 {} && ".format( output.unmapped, output.unmapped_r1, output.unmapped_r2 )

        # Pair reads if they are amplicons, so paired information can inform contamination.
        if itype == "amplicon":
            merged = "{}/_barcode/{}.assembled.fastq".format( out_dir, wilcards.sample )
            merged_prefix = "{}/_barcode/{}".format( out_dir, wilcards.sample )
            command += "pear -n 356 -m 398 -f {output.unmapped_r1} -r {output.unmapped_r2} -o {out_dir}/_barcode/{wildcards.sample} && ".format( output.unmapped_r1, output.unmapped_r2, merged_prefix )
            alignment_input = [ output.merged ]
            alignment_reference = config["barcode_check"]["barcodes"]
        else:
            alignment_input = [output.unmapped_r1, output.unmapped_r2]
            alignment_reference = config["barcode_check"]["barcode_regions"]

        # Align the unmapped reads to to barcode reference. This differs based on whether the reads where generated
        # using nextera or an amplicon based method. Previous versions used an alignment score limit of 346 for amplicon
        command += "bwa mem -T 346 {} {} | samtools view -F 4 -Sbu | samtools sort -o {} && ".format( alignment_reference, " ".join( alignment_input ), output.barcode_bam )
        command += "samtools index {} && ".format( output.barcode_bam )
        command += "python3 {config[scripts]}/contamination.py -l {} -b {} {} {} {}".format( limit, alignment_reference, "-n" if itype=="nextera" else "-a", output.barcode_bam, output.barcode_stats )
        call( command, shell=True )
        
rule align_reads:
    input:
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][0],
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][1]
    output:
        "{out_dir}/_untrimmed/{sample}.aligned.sorted.bam"
    shell:
        "bwa mem {ref_sequence} {input[0]} {input[1]} | samtools view -F 4 -Sbu | samtools sort -o {output} &&"
        "samtools index {output}"
