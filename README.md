# locuszoomplot

Script to create locuszoomplot with GWAS and fine-mapping results.

## Install all needed R packages:
Rscript scripts/Install_packages.R 

## Example usage for showing the help page:
Rscript scripts/locuszoomplot.R -h


Options available:


Mandatory fields:

--file=
		Results files from susie (example_files/FEMALE_INFERT_NEW.chr6.149730382-152730382.susie.snp).

--annofile=
		Variant annotation file including most severe consequence and rsid columns .

--geneFile=
		File with gene information (example_files/Gene_data_for_locuszoomplot.txt).

--PlotName=
		Prefix for the plots name.


Optional fields:

--Nrows=
		Number of rows for genes. Default is 4.

--PlotTitle=
		Title to give for a plot (will be plotted in italics). Default is empty.

--proteinCodingOnly=
		Whether to print only protein coding genes. Default is TRUE.

-h, --help
		Show this help message and exit
	

	
