# locuszoomplot

Script to create locuszoomplots with GWAS and fine-mapping results.

Example usage for the help page:
Rscript scripts/locuszoomplot.R -h


Options available:


Mandatory fields:

--file=FILE
		Results files from susie (example_files/FEMALE_INFERT_NEW.chr6.149730382-152730382.susie.snp).

--annofile=ANNOFILE
		Variant annotation file including most severe consequence and rsid columns .

--geneFile=GENEFILE
		File with gene information (example_files/Gene_data_for_locuszoomplot.txt).

--PlotName=PLOTNAME
		Prefix for the plots name.


Optional fields:

--Nrows=NROWS
		Number of rows for genes. Default is 4.

--PlotTitle=PLOTTITLE
		Title to give for a plot (will be plotted in italics). Default is empty.

--proteinCodingOnly=PROTEINCODINGONLY
		Whether to print only protein coding genes. Default is TRUE.

-h, --help
		Show this help message and exit
	

	
