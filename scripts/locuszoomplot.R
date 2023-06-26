#!/usr/bin/env Rscript
options(stringsAsFactors=F)

## load R libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(httr)
library(xml2)
library(rlist)
library(jsonlite)
options(scipen = 999) 

option_list <- list(
  make_option("--file", type="character",default="",
              help="Results files from susie"),
  make_option("--annofile", type="character",default="",
              help="Variant annotation file including most severe consequence and rsid."),
  make_option("--geneFile", type = "character", default="",
              help="File with gene information."),
  make_option("--PlotName", type = "character", default="",
              help="Prefix for the plots name."),
  make_option("--Nrows", type = "integer", default="4",
              help="Number of rows for genes. Optional. Default is 4."),
  make_option("--PlotTitle", type = "character", default = "",
              help = "Title to give for a plot (will be plotted in italics). Optional. Default is empty."),
  make_option("--proteinCodingOnly", type = "logical", default="T",
              help="Whether to print only protein coding genes. Optional. Default is TRUE."))

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

break_points <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors <- c('#01007B',
            '#96CCF5',
            '#75FB4C',
            '#F2A83B',
            '#EA3423')

d<-fread(opt$file)
chrom= gsub("chr", "", d$chromosome[1])
max<-max(-log10(d$p))+0.2*max(-log10(d$p))
lead_index=which(d$p==min(d$p))
d$variation2<-gsub("_", ":", d$rsid)
d$variation2<-gsub("chr", "", d$variation2)
var<-d$variation2[lead_index]

minbp<-min(d$position)-10000
maxbp<-max(d$position)+10000
win<-maxbp-minbp

cmd = paste0("resp<-GET('http://api.finngen.fi/api/ld?variant=",var,"&window=",win,"&panel=sisu42&r2_thresh=0')")
eval(parse(text=cmd))
jsonRespParsed<-content(resp, as="parsed") 
modJson<-jsonRespParsed$ld
lddat<-list.stack(modJson)
lddat$ld<-lddat$r2
lddat<-select(lddat, c("variation2", "ld"))

d<-left_join(d, lddat)
d$ld[lead_index]<-1

a<-fread(opt$annofile)
a<-a[a$variant%in%d$rsid,]
names(d)[which(names(d)=="rsid")]<-"variant"
d<-left_join(d, a)

coding<-c("transcript_ablation", "splice_donor_variant", "stop_gained", "splice_acceptor_variant", "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant")
d$col = colors[ as.numeric( cut( d$ld, breaks = break_points ) ) ]
d$col<-ifelse(d$ld==0, '#01007B', d$col)
d$col[lead_index]<-'#722CC5'

d$shape = ifelse(d$most_severe%in%coding, 22, 21)
d$shape[lead_index] = 23
d$size <- 2
d$size[lead_index] <- 2.8

n_cs<-length(table(d$cs))-1
if (n_cs>=4) {
  print("More than 4 credible sets, credible sets beyond 4th won't be displayed.")
}
pchs<-c(23, 22, 24, 25)
pchs_to_legend= c(pchs[1:n_cs], 21)
legends<-c()
for(i in 1:n_cs){
  legends[i]<-paste0("In credible set ",i,"")
}
legends[n_cs+1] <-"Not in any credible set"

d$shape_cs<-ifelse(d$cs==1, 23, ifelse(d$cs==2, 22, ifelse(d$cs==3, 24, ifelse(d$cs==4, 25, 21))))

max_prob = max(d$prob)+0.1*max(d$prob)
rsid_lead=d$rsid[lead_index]

gene<-fread(opt$geneFile)

if(isTRUE(opt$proteinCodingOnly)){
  gene<-gene[gene$Gene_type=="protein_coding",]
}

minbp<-min(d$position)
maxbp<-max(d$position)
gene<-gene[gene$Chromosome==chrom & gene$start>minbp & gene$end<maxbp,]
gene<-subset(gene, gene$Gene_name!="")
gene<-gene[order(gene$start),]
gene$namesplace<-(gene$start+gene$end)/2

n_genes<-nrow(gene)
n_rows=opt$Nrows
gene$row <- c(rep(1:n_rows, length.out = n_genes))
gene$direction<-gene$Strand
GENE_plot<-(paste0("'", gene$Gene_name, "'"))

gene$gene_name<-ifelse(gene$direction=="-1",
                       paste0("expression('' %<-% italic(", GENE_plot, "))"),
                       paste0("expression(italic(", GENE_plot, ") %->% '')"))

cmd = paste0("tiff('",opt$PlotName,"_",rsid_lead,".tiff', 14, 12, res = 300, units = 'in', bg = 'white')")
eval(parse(text = cmd))
par(fig=c(0,10,5,10)/10)
par(mar=c(0.5,5,4,5))
cmd = paste0("plot(d$position/1000000, -log10(d$p),
     ylim = c(0, max),
     bg = d$col,
     xaxs='i',
     yaxs='i',
     las = 1,
     main = substitute(paste(italic('",opt$PlotTitle,"'))),
     cex.main = 1.5,
     pch = d$shape,
     cex = d$size,
     xaxt='n',
     ylab = '-log10(P-value)',
     xlab = '')")
eval(parse(text = cmd))
abline(h = -log10(5E-8), lty = 2, col = "grey40")
points(d$position[d$col=='#722CC5']/1000000, -log10(d$p[d$col=='#722CC5']),
       bg = d$col[d$col=='#722CC5'],
       pch = d$shape[d$col=='#722CC5'],
       cex = d$size[d$col=='#722CC5'])
points(d$position[d$shape==22]/1000000, -log10(d$p[d$shape==22]),
       bg = d$col[d$shape==22],
       pch = d$shape[d$shape==22],
       cex = d$size[d$shape==22])
cmd = paste0("text(d$position[lead_index]/1000000,
     -log10(d$p[lead_index])+0.5,
     '",rsid_lead,"',
     pos = 3,
     cex = 1.2)")
eval(parse(text=cmd))
cmd = paste0("legend( 'topleft',
        inset=.01,
        legend = paste0( break_points[1:(length(break_points) -1) ], ' - ',
                         break_points[2:length(break_points) ] ),
        fill = colors,
        cex = 1.5,
        pt.cex = 3,
        title = expression(paste('r'^'2','(",rsid_lead,")')))")
eval(parse(text=cmd))
legend('topright',
       pch = c(19, 15),
       cex = 1.7,
       legend = c("Non-coding", "Coding"))

par(new = T)
par(fig=c(0,10,2,5)/10)
par(mar=c(5.1,5,0,5))
plot(d$position/1000000, d$prob,
     bg = d$col,
     ylim = c(0, max_prob),
     xaxs="i",
     yaxs="i",
     las = 1,
     pch = d$shape_cs,
     cex = d$size*0.7,
     ylab = "Posterior Incusion Probability (PIP)",
     xlab =paste0("Position in chromosome ",chrom," (Mb)"))
legend('topleft',
       pch = pchs_to_legend,
       legend = legends,
       cex = 1.5,
       pt.bg = "darkblue")
par(new = T)
par(fig=c(0,10,0,2)/10)
par(mar=c(0.5,5,0.5,5))
plot(gene$start/10000000, n_rows+1-gene$row,
     xlim=c(minbp/10000000, maxbp/10000000),
     ylim = c(0.8, n_rows+0.4),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     type="n")
segments(gene$start/10000000, n_rows+1-gene$row,
         gene$end/10000000, n_rows+1-gene$row,
         col = "#01007B", lty = 1, lwd = 2,
         xlim=c(minbp/10000000, maxbp/10000000))
for(i in 1:nrow(gene)){
  text(gene$namesplace[i]/10000000, n_rows+1-gene$row[i],
       labels = eval(parse(text = gene$gene_name[i])),
       pos = 3, font = 3,
       col =  "#01007B", cex = 0.7)
}
dev.off()
