import itertools

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import log10

def map_linewidth( value, upper, lower, upper_prime, lower_prime ):
    r = ( upper_prime - lower_prime ) / ( upper - lower )
    return ( value - lower ) * r + lower_prime

def plot_contamination( contamination_events, output, exclusion ):
    plt.figure( dpi=200, figsize=(10,5) )
    ax = plt.subplot()

    # Draw lines for plate.
    for i in range( 13 ) :
        if i == 0 :
            ls = "solid"
        else :
            ls = "dashed"

        if i == 12 : ls = "solid"
        ax.vlines( i + 0.5, 0.5, 8.5, zorder=1, linestyles=ls )

        if i < 9 :
            if i == 8 :
                ls = "solid"
            ax.hlines( i + 0.5, 0.5, 12.5, zorder=1, linestyles=ls )

    # Draw the contamination events to the image.
    arrowstyle = "Simple,tail_width=0.5,head_width=4,head_length=8"
    kw = dict( arrowstyle=arrowstyle, color="k" )
    for i in contamination_events:
        ax.add_patch( patches.FancyArrowPatch( i[0], i[1], linewidth=map_linewidth( log10( i[2] ), 1.69, -4, 4, 0 ), connectionstyle="arc3,rad=.1", **kw,  ) )

    # Draw the excluded boxes to the image.
    for i in exclusion:
        new_i = (i[0]-0.5, i[1]-0.5)
        ax.add_patch( patches.Rectangle( new_i, 1, 1, color="#BEF2E6", zorder=0 ))

    # Set up y ticks
    ax.set_yticks( range( 1, 9 ) )
    ax.set_yticklabels( ["H", "G", "F", "E", "D", "C", "B", "A"] )
    ax.set_ylim((0.49, 8.51 ))

    # Set up xticks
    ax.set_xticks( range( 1, 13 ) )
    ax.xaxis.tick_top()
    ax.set_xlim( (0.49, 12.51) )

    ax.spines['right'].set_visible( False )
    ax.spines['top'].set_visible( False )
    ax.spines['bottom'].set_visible( False )
    ax.spines['left'].set_visible( False )

    plt.tight_layout()
    plt.savefig( output )

def guess_contamination_events( matrix, inclusion_limit, contamination_limit ):

    # Generate barcode x sample dictionary
    barcode_dict = matrix.idxmax().to_dict()

    # return object; format is contaminant_source, contaminant_sink, magnitude
    events = list()

    # generate proportional matrix
    matrix_pro = matrix.apply( lambda x: x / x.sum() )

    for i in matrix_pro.loc[:,matrix.sum() > inclusion_limit]:
        temp = matrix_pro.loc[matrix_pro[i] > contamination_limit,i]
        temp = temp.drop( barcode_dict[i] )
        #print( "{}:{}".format( i, len( temp ) ) )
        for j in temp.iteritems():
            events.append( ( get_position(j[0]), get_position(barcode_dict[i]), j[1], j[0], barcode_dict[i] ) )

    exclusion_list = [get_position(barcode_dict[i]) for i in matrix.columns[matrix.sum() < inclusion_limit]]
    return events, exclusion_list

def get_majority( series ):
    forward = series.loc[series["Forward"],"Barcode"].iloc[0]
    reverse = series.loc[~series["Forward"],"Barcode"].iloc[0]
    return forward, reverse

def get_all_barcodes( series ):
    permutations = itertools.product( series.loc[series["Forward"],"Barcode"], series.loc[~series["Forward"],"Barcode"] )
    return [i for i in permutations]

def get_position( bc_name, columnwise=True ):
    bc_number = int( bc_name.split( "-" )[-1] )
    x_pos = (bc_number + 1) % 12
    if x_pos == 0:
        x_pos = 12
    y_pos = 8 - int( bc_number / 12 )

    return x_pos, y_pos

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Takes a sample x barcode read matrix and estimates graphs estimated contamination events." )

    # Initialize positional arguments
    parser.add_argument( "matrix", help="sample by barcode read matrix" )
    parser.add_argument( "output", help="Ouput image save location" )

    # Intialized optional arguments
    order_group = parser.add_mutually_exclusive_group( required=False )
    order_group.add_argument( "-c", "--columnwise", help="order of samples in matrix or orderfile is columnwise (i.e A1, A2, A3, etc)", action="store_true" )
    order_group.add_argument( "-r", "--rowwise", help="order of samples in matrix or orderfile is rowwise (i.e. A1, B1, C1, etc)", action="store_true" )

    parser.add_argument( "-l", "--limit", default=10000, type=float, help="limits the inclusive of sample in output graph. Values above one represent a read limit, while values below one represent presentile limit" )
    parser.add_argument( "-L", "--contaminant-limit", default=0.01, type=float, help="Limit above which contaminant event is included in output graph. Values above one represent a read limit, while values below one represent a percentage limit" )

    args = parser.parse_args()

    # This is temporary, and a sorted matrix should be used instead.
    bc = pd.read_csv( args.matrix, index_col="Barcode" )
    bc = bc.reindex( sorted( bc.columns ), axis=1 )
    bc = bc.drop( columns=["CSUNegCtrl", "Undetermined"] )

    ce,el = guess_contamination_events( bc, args.limit, args.contaminant_limit )


    plot_contamination( ce, args.output, el )

    for i in ce: print(i)

    #barcode_dict = dict()
    #for l, i in enumerate( range( 1, 9 ) ) :
    #    for k, j in enumerate( range( 9, 21 ) ) :
    #        pairing = ("Barcode_Forward_{}".format( i ), "Barcode_Reverse_{}".format( j ))
    #        barcode_dict[pairing] = "BarcodeOligo-{}".format( (i - 1) * 12 + k )
#
    #data = { "Barcode" : ["Barcode_Forward_1", "Barcode_Forward_2", "Barcode_Forward_3", "Barcode_Reverse_4", "Barcode_Reverse_5"],
    #         "Forward" : [True, True, True, False, False],
    #         "A" : [50, 5, 0, 50, 5],
    #         "B" : [50, 5, 0, 0, 55],
    #         "C" : [0, 50, 0, 50, 0],
    #         "D" : [50, 0, 0, 50, 0],
    #         "E" : [0, 0, 50, 50, 0] }
    #data = pd.DataFrame( data )
    #data = data.set_index( ["Barcode", "Forward"] )
    #print( data.head() )
#
    #data["Forward","A"]