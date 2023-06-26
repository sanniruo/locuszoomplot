# locuszoomplot

Script to create locuszoomplot with GWAS and fine-mapping results.

### Install all needed R packages:
`Rscript scripts/Install_packages.R`

### Example usage for showing the help page:
`Rscript scripts/locuszoomplot.R -h`


## Parameters available:

### Mandatory arguments

```
  --file=FILE
    Results files from susie (example_files/FEMALE_INFERT_NEW.chr6.149730382-152730382.susie.snp).

  --annofile=ANNOFILE
    Variant annotation file including most severe consequence and rsid.

  --geneFile=GENEFILE
    File with gene information (example_files/Gene_data_for_locuszoomplot.txt).

  --PlotName=PLOTNAME
    Prefix for the plots name.

```

### pptional arguments:

```
  -h, --help            show this help message and exit
    --Nrows=NROWS
    Number of rows for genes. Optional. Default is 4.

  --PlotTitle=PLOTTITLE
    Title to give for a plot (will be plotted in italics). Optional. Default is empty.

  --proteinCodingOnly=PROTEINCODINGONLY
    Whether to print only protein coding genes. Optional. Default is TRUE.
```
	

	
