startTime=Sys.time()
#!/usr/bin/env Rscript
options(stringsAsFactors=F)

## load R libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(XLConnect)
library(httr)
library(xml2)
library(rlist)
library(jsonlite)
options(scipen = 999) 

option_list <- list(
  make_option("--file", type="character",default="",
              help="Results files from susie"),
  make_option("--annofile", type="character",default="",
              help="Variant annotation file"),
  make_option("--gene_data", type = "character", default="",
              help="File with gene information, downloaded from UCSC"))

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
d$shape_cs<-ifelse(d$cs==1, 23, ifelse(d$cs==2, 22, 21))
max_prob = max(d$prob)+0.2*max(d$prob)
rsid_lead=d$rsid[lead_index]

gene<-fread(opt$gene_data)
GENE<-unique(gene$hg38.kgXref.geneSymbol)
GENE<-GENE[!grepl("ENSG0", GENE)]
gene_data<-data.frame(GENE=GENE)

for(i in 1:length(GENE)){
  sub<-subset(gene,  gene$hg38.kgXref.geneSymbol==GENE[i])
  gene_data$start[i] <- min(sub$hg38.knownGene.txStart)
  gene_data$end[i] <-max(sub$hg38.knownGene.txEnd)
  gene_data$direction[i]<-sub$hg38.knownGene.strand[1]
}
gene_data$namesplace<-(gene_data$start+gene_data$end)/2

gene_data$row2 <- c(rep(1:5, times = 9))
GENE_plot<-(paste0("'", GENE, "'"))

gene_data<-gene_data[order(gene_data$start),]
gene_data$gene_name<-ifelse(gene_data$direction=="-",
                            paste0("expression('' %<-% italic(", GENE_plot, "))"),
                            paste0("expression(italic(", GENE_plot, ") %->% '')"))

genenames<-c()
starts<-c()
ends<-c()
for(i in 1:length(GENE)){
  name<-GENE[i]
  sub<-gene[gene$hg38.kgXref.geneSymbol==name,]
  n<-length(unlist(strsplit(sub$hg38.knownGene.exonStarts, ",")))
  genename <- rep(name, times = n)
  start <- unlist(strsplit(sub$hg38.knownGene.exonStarts, ","))
  end <- unlist(strsplit(sub$hg38.knownGene.exonEnds, ","))
  genenames<-c(genenames, genename)
  starts<-c(starts, start)
  ends<-c(ends, end)
}

tick_data<-data.frame(GENE = genenames, start = starts, end = ends)
m<-match(tick_data$GENE, gene_data$GENE)
tick_data$row<-gene_data$row2[m]
tick_data$start<-as.numeric(tick_data$start)
tick_data$end<-as.numeric(tick_data$end)


pdf("SYNE1_testi.pdf", 14, 12, bg = "white")
par(fig=c(0,10,5,10)/10)
par(mar=c(0.5,5,4,5))
plot(d$position/1000000, -log10(d$p),
     ylim = c(0, max),
     bg = d$col,
     xaxs="i",
     yaxs="i",
     las = 1,
     pch = d$shape,
     cex = d$size,
     xaxt='n',
     ylab = "-log10(P-value)",
     xlab = "")
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
       pch = c(21, 22),
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
legend('topright',
       pch = c(23, 22, 21),
       legend = c("In CS #1", "In CS #2", "Not in any CS"))
par(new = T)
par(fig=c(0,10,0,2)/10)
par(mar=c(0.5,5,0.5,5))
plot(gene_data$start/10000000, 6-gene_data$row2,
     xlim=c(149730382/10000000,152730382/10000000),
     ylim = c(0.8, 5.4),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     type="n")
segments(gene_data$start/10000000, 6-gene_data$row2,
         gene_data$end/10000000, 6-gene_data$row2,
         col = "#01007B", lty = 1, lwd = 2,
         xlim=c(149730382/10000000,152730382/10000000))
for(i in 1:nrow(gene_data)){
  text(gene_data$namesplace[i]/10000000, 6-gene_data$row2[i],
       labels = eval(parse(text = gene_data$gene_name[i])),
       pos = 3, font = 3,
       col =  "#01007B", cex = 0.5)
}
segments(tick_data$start/10000000, 6-tick_data$row-0.02,
         tick_data$start/10000000, 6-tick_data$row+0.02,
         col = "#01007B", lty = 1, lwd = 2,
         xlim=c(149730382/10000000,152730382/10000000))
dev.off()

endTime=Sys.time()
endTime-startTime