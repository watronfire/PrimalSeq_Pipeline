import os
from fnmatch import fnmatch
import pandas as pd
from subprocess import call, check_output
from Bio import SeqIO
import itertools 

# Inputs and required files for pipeline
in_dir=config["in_dir"]
out_dir=config["out_dir"]
ref_sequence=config["ref_sequence"]
barcodes_f=config["barcode_check"]["forward_barcode"]
barcodes_r=config["barcode_check"]["reverse_barcode"]
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

# Helper function to identify barcodes from regions
def extract_barcode( barcode_string ):
    parts = barcode_string.split( "-" )
    forward = int( parts[0][2:] )
    reverse = int( parts[1][2:] )
    
    return( "BarcodeOligo-{}".format( (forward - 1) * 12 + (reverse - 8) ) )


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
        #expand( "{out_dir}/_consensus/{sample}.fa", out_dir = out_dir, sample = SAMPLES ),
        #expand( "{out_dir}/_barcode/{sample}_bs.bam", out_dir=out_dir, sample=SAMPLES ),
        #expand( "{out_dir}/_final/{full_name}.fasta", out_dir=out_dir, full_name=final_file.values() ),
        #expand( "{out_dir}/_coverage/{sample}.coverage.png", out_dir = out_dir, sample = SAMPLES ),
        expand( "{out_dir}/_untrimmed/{sample}.aligned.sorted.bam", out_dir = out_dir, sample = SAMPLES ),
        expand( "{out_dir}/_barcode/{sample}.csv", out_dir=out_dir, sample=SAMPLES ),
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

        bc_dfs = list()

        # Alignment Statistics
        data_dict = { "Sample" : list(), "Total_Reads" : list(), "Percent_Coverage" : list(), "Mean_Coverage" : list(),"Barcode_Reads" : list(),  "Barcode" : list(),  "Contamination" : list() }
        for al, bc  in zip( sorted( input.align ), sorted( input.barcodes ) ):
            sample = os.path.splitext( os.path.basename( al ) )[0]
            data_dict["Sample"].append( sample )
            data_dict["Total_Reads"].append( int( check_output( ["samtools", "view", "-c", al] ).strip() ) )
            
            coverage = check_output( ["samtools", "depth", "-aa", "-d", "0", al] )

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
            
            # Get percent contamination
            bc_df = pd.read_csv( bc )
            bc_df = bc_df.set_index( "Barcode" )
            bc_dfs.append( bc_df )
            data_dict["Barcode"].append(  bc_df.idxmax().item() )
            data_dict["Barcode_Reads"].append( bc_df.sum().item() )
            data_dict["Contamination"].append( 1.0 - ( bc_df.max() / bc_df.sum() ).item() )
            
        data_df = pd.DataFrame( data_dict )
        data_df.to_csv( str( output.alignment_stats ) )

        # Barcode Stats
        bc_all = pd.concat( bc_dfs, axis=1 )
        bc_all.to_csv( output.barcode_stats )
        
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

if config["barcode_check"]["type"] == "nextera":
    rule barcode_check_nextera:
        input:
            unmapped = "{out_dir}/_untrimmed/{sample}.unmapped.bam"
        output:
            unmapped_r1 = temp( "{out_dir}/_barcode/{sample}_r1.fastq" ),
            unmapped_r2 = temp( "{out_dir}/_barcode/{sample}_r2.fastq" ),
            barcode_bam = "{out_dir}/_barcode/{sample}_bc.bam",
            completed = "{out_dir}/_barcode/{sample}.csv"
        run:
            command = "module load bedtools && samtools index {} &&".format( input.unmapped )
            command += "bedtools bamtofastq -i {} -fq {} -fq2 {} && ".format( input.unmapped, output.unmapped_r1, output.unmapped_r2 )
            call( command, shell=True )

            # Align reads to barcode regions
            command = "bwa mem {} {} {} | samtools view -F 4 -Sbu | samtools sort -o {} && ".format( config["barcode_check"]["barcode_regions"], output.unmapped_r1, output.unmapped_r2, output.barcode_bam )
            command += "samtools index {} && ".format( output.barcode_bam )
            command += "python3 {config[scripts]}/contamination.py -l {} -b {} -n {} {}".format( limit, config["barcode_check"]["barcode_regions"], output.barcode_bam, output.completed )
            call( command, shell=True )

else:
    rule barcode_check_amplicon:
        input:
            unmapped = "{out_dir}/_untrimmed/{sample}.unmapped.bam"
        output:
            unmapped_r1 = temp( "{out_dir}/_barcode/{sample}_r1.fastq" ),
            unmapped_r2 = temp( "{out_dir}/_barcode/{sample}_r2.fastq" ),
            completed = "{out_dir}/_barcode/{sample}.csv"
        run:
            # Load required modules
            command = "module load bedtools && "
            command += "samtools index {} && ".format( input.unmapped )

            # Seperate read pairs into seperate files from alignment
            command += "bedtools bamtofastq -i {} -fq {} -fq2 {}".format( input.unmapped, output.unmapped_r1, output.unmapped_r2 )

            call( command, shell=True )

            # Demultiplex with cutadapt assuming regular orientationi
            round1_1 = os.path.join( out_dir, "_barcode/r1.{}.{{name1}}-{{name2}}.1.fastq.gz".format( wildcards.sample ) )
            round1_2 = os.path.join( out_dir, "_barcode/r1.{}.{{name1}}-{{name2}}.2.fastq.gz".format( wildcards.sample ) )
            command = "cutadapt --quiet --action=none -g file:{} -G file:{} -o {} -p {} {} {}".format( barcodes_f, barcodes_r, round1_1, round1_2, output.unmapped_r1, output.unmapped_r2 )
            call( command, shell=True )   

            unknown_r1 = os.path.join( out_dir, "_barcode/r1.{}.unknown-unknown.1.fastq.gz".format( wildcards.sample ) )
            unknown_r2 = os.path.join( out_dir, "_barcode/r1.{}.unknown-unknown.2.fastq.gz".format( wildcards.sample ) )

            # Demultiplex with cutadapot assuming goofy orientation.
            round2_1 = os.path.join( out_dir, "_barcode/r2.{}.{{name1}}-{{name2}}.1.fastq.gz".format( wildcards.sample ) )
            round2_2 = os.path.join( out_dir, "_barcode/r2.{}.{{name1}}-{{name2}}.2.fastq.gz".format( wildcards.sample ) )
            command = "cutadapt --quiet --action=none -g file:{} -G file:{} -o {} -p {} {} {} &&".format( barcodes_f, barcodes_r, round2_2, round2_1, unknown_r2, unknown_r1 )
            call( command, shell=True )

            # Identify output files of cutadapt
            files = [i for i in os.listdir( os.path.join( out_dir, "_barcode" )  ) if fnmatch( i, "*.{}.*.1.fastq.gz".format( wildcards.sample ) ) and "unknown" not in i]
            
            # lambda function to sort files generated by cutadapt
            combiner = lambda x : x.split( "." )[2]
            files = sorted( files, key=combiner )

            # group files by barcode 
            with open( output.completed, "w" ) as completed:
                completed.write( "Barcode,{}\n".format( wildcards.sample ) )
                for k, g in itertools.groupby( files, key=combiner ):
                    file_loc = [os.path.join( out_dir, "_barcode/{}".format( i ) ) for i in g]
                    _ = check_output( "echo $(zcat {} | wc -l )/4 | bc ".format( " ".join( file_loc ) ), shell=True )
                    bc = extract_barcode( k )
                    completed.write( "{},{}\n".format( bc, int( _ ) ) )
            
            # Delete demultiplexed read files.
            call( "rm {}/_barcode/*.{}.*.fastq.gz".format( out_dir, wildcards.sample ), shell=True )  

rule align_reads:
    input:
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][0],
        lambda wildcards: "{in_dir}".format( in_dir = in_dir ) + "/" + SAMPLES[wildcards.sample][1]
    output:
        temp( "{out_dir}/_untrimmed/{sample}.alignedOnly.bam" ),
        "{out_dir}/_untrimmed/{sample}.unmapped.bam",
        "{out_dir}/_untrimmed/{sample}.aligned.sorted.bam"
    shell:
        "bwa mem {ref_sequence} {input[0]} {input[1]} > {output[0]} &&"
        "samtools view -F 4 -Sbu {output[0]} | samtools sort -o {output[2]} &&"
        "samtools index {output[2]} &&"
        "samtools view -b -f 12 -F 256 {output[0]} | samtools sort -o {output[1]}"
