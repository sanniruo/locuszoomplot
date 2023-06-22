# locuszoomplot

Script to create locuszoomplots with GWAS and fine-mapping results.

Example usage for the help page:
Rscript scripts/locuszoomplot.R -h

Options available:

Mandatory fields:

--files	Results files from SUSIE. (example_files/FEMALE_INFERT_NEW.chr6.149730382-152730382.susie.snp)

--annofile	Variant annotation file including most severe consequence and rsid.

--geneFile	File with gene information. (example_files/Gene_data_for_locuszoomplot)

--PlotName	Prefix for the plot's name.

Optional fields:

--Nrows	The number of rows for genes. Default is 4. (optional)

--proteinCodingOnly	Whether to print only protein-coding genes. The default is TRUE. (optional)

--help Show this help message and exit
