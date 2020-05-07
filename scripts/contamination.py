import argparse
import pysam
import pandas as pd
import os

def index_stat( alignment, barcodes, limit ):
    for i in alignment.fetch():
        if i.get_tag( "AS" ) > limit:
            barcodes.loc[i.reference_name, "Count" ] += 1
    return barcode_df

def check_input( file ):
    loaded_file = pysam.AlignmentFile( file )
    assert loaded_file.has_index()
    return loaded_file

def generate_barcode_df( file ):
    with open( file ) as bc_fasta:
        barcodes = [i[1:].strip() for i in bc_fasta if ">" in i]
    return pd.DataFrame( { "Barcode" : barcodes,
                           "Count" : [0] * len( barcodes ) } )

if __name__ == "__main__":

    # Initialize argument parser
    parser = argparse.ArgumentParser( description="Takes a barcode-aligned bam and returns contamination statistics" )

    # Initialize positional arguments
    parser.add_argument( "alignment", help="input alignment, containing reads mapped to barcodes" )
    parser.add_argument( "output", help="output csv" )

    # Initialize optional arguments
    parser.add_argument( "-l", "--limit", help="minimum required alignment score to filter accurate reads", type=int, required=True )
    parser.add_argument( "-b", "--barcodes", help="fasta of barcode files" )

    group = parser.add_mutually_exclusive_group( required=True )
    group.add_argument( "-n", "--nextera", help="specify nextera-based reads", action="store_true" )
    group.add_argument( "-a", "--amplicon", help="specify amplicon reads", action="store_true" )

    # Parse arguments
    args = parser.parse_args()

    # Check alignment file
    try:
        loaded_alignment = check_input( args.alignment )
    except AssertionError:
        print( "Error: Input alignment must be indexed. Call `samtools index <alignment>` first, next time." )
        exit( 42 )

    # Generate barcode_dict
    barcode_df = generate_barcode_df( args.barcodes )
    if args.nextera:
        barcode_df["Forward"] = True
        barcode_df.loc[barcode_df["Barcode"].str.contains( "Reverse" ), "Forward"] = False

    barcode_df = barcode_df.set_index( "Barcode" )

    # Add idxstats to barcode dataframe.
    barcode_df = index_stat( loaded_alignment, barcode_df, args.limit )

    # Rename count column so merge can occur later
    name = os.path.basename( args.alignment )
    name = str( os.path.splitext( name )[0] )
    barcode_df = barcode_df.rename( columns={"Count" : name } )

    barcode_df.to_csv( args.output )
