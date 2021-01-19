import matplotlib
#matplotlib.use( "Agg" )
import matplotlib.pyplot as plt
import argparse
import pandas as pd


def plot_coverage( coverage, coverage_average, input_name, output, bed=None ):
	plt.figure( dpi=200, figsize=(10,4) )
	ax = plt.subplot()

	ax.plot( "Position", "Coverage", data=coverage.reset_index(), color="#49BC84", zorder=3 )
	ax.fill( "Position", "Coverage", data=coverage.reset_index(), color="#A4DDC1", zorder=2 )
	ax.set_yscale( "log" )

	ax.axhline( coverage_average, 0, 1, color="black", linestyle="dashed", zorder=4 )
	ax.grid( True, axis="y", linewidth=0.5, color="#969696", zorder=1 )

	if bed is not None:
		for name, row in bed.iterrows():
			temp = coverage[row.start:row.end]
			if any( temp["Coverage"] < 15 ):
				ax.axvspan( row.start, row.end, zorder=0, color="#dfdfdf" )
				ax.text( ( row.start + row.end ) / 2, coverage["Coverage"].max(), name + 1,
						 horizontalalignment="center",
						 fontdict={"fontsize" : "x-small", "fontweight" : "bold"} )


	ax.set_ylabel( "Reads", weight="bold" )
	ax.set_xlabel( "Genome Position (bp)", weight="bold" )
	ax.set_title( "Coverage Map for " + input_name, weight="bold" )

	plt.tight_layout()
	plt.savefig( output )


def graph_coverage( args ):
	input_name = args.input.split( "/" )[-1].split( "." )[0]
	coverage = pd.read_csv( args.input, header=None, names=["Reference", "Position", "Coverage"], sep="\t" )
	coverage = coverage.set_index( "Position" )
	coverage_average = coverage["Coverage"].mean()

	#print( coverage_average )
	#print( coverage.head() )

	inserts = pd.read_csv( args.bed, sep="\t", header=None, names=["Reference", "start", "end", "name", "pool", "strand"] )

	plot_coverage( coverage, coverage_average, input_name, args.output, bed=inserts )

if __name__ == "__main__":
	parser = argparse.ArgumentParser( description="generates a plot showing the coverage per nucleotide across the genome" )

	# Initialize positional arguments
	parser.add_argument( "-i", "--input", required=True, help="TSV file containing the coverage per position" )
	parser.add_argument( "-b", "--bed", help="bed file containing the position of primers" )
	parser.add_argument( "-o", "--output", required=True, help="location to save image" )

	arguments = parser.parse_args()

	graph_coverage( arguments)
