# locuszoomplot

Script to create locuszoomplots with GWAS and fine-mapping results.

Rscript scripts/locuszoomplot.R -h

Options available:


	--file=FILE
		Results files from susie.

	--annofile=ANNOFILE
		Variant annotation file including most severe consequence and rsid.

	--geneFile=GENEFILE
		File with gene information.

	--PlotName=PLOTNAME
		Prefix for the plots name.

	--Nrows=NROWS
		Number of rows for genes. Default is 4. (optional)

	--proteinCodingOnly=PROTEINCODINGONLY
		Whether to print only protein coding genes. Default is TRUE. (optional)

	-h, --help
		Show this help message and exit
