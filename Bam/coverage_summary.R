################################################################################
#Reports some statistics on a coverage file
#
#Depth file created with:
#   samtools depth -a <input.bam>
#
#@author:charles.hefer@agresearch.co.nz
################################################################################
nz_repo = c("https://cran.stat.auckland.ac.nz/")

if (!require("tidyverse")) {
  install.packages("tidyverse", repos = nz_repo)
}
library("tidyverse")

if (!require("cowplot")) {
  install.packages("cosplot", repos = nz_repo)
}
library("cowplot")



args = commandArgs(trailingOnly=TRUE)

#for debuggins
args <- c("20-PoM16-1_vs_DUG42.depth", 
          "01-ViL-Nal1_vs_DUG42.png",
          "01-ViL-Nal1_vs_DUG42.summary")

if (length(args)<3) {
  stop("Specify both the input depth file, output png and output summary files")
}

threshold <- function(df) {
  th_df <- NULL
  th_df$x0 <- sum(df$V3==0)/length(df$V3)
  th_df$x1 <- sum(df$V3==1)/length(df$V3)
  th_df$x2 <- sum(df$V3==2)/length(df$V3)
  th_df$x3 <- sum(df$V3==3)/length(df$V3)
  th_df$x4 <- sum(df$V3==4)/length(df$V3)
  th_df$x5 <- sum(df$V3==5)/length(df$V3)
  th_df$x6to10 <- sum(df$V3>5 & df$V3<=10)/length(df$V3)
  th_df$x11to20 <- sum(df$V3>10 & df$V3<=20)/length(df$V3)
  th_df$x21to50 <- sum(df$V3>20 & df$V3<=50)/length(df$V3)
  th_df$x51to100 <- sum(df$V3>50 & df$V3<=100)/length(df$V3)
  th_df$gt100 <- sum(df$V3>100)/length(df$V3)
  return(data.frame(th_df))
}

fh <- read.table(args[1])
p1 <- fh %>%
  ggplot(aes(x=V3)) +
  geom_histogram(aes(y=after_stat(count/sum(count))), binwidth=50) + 
  scale_y_continuous(labels=scales::percent) + 
  theme_bw() + 
  ylab("Percentage of genome covered") + 
  xlab("Depth of coverage (max 8000x)") + 
  ggtitle(args[1])

p2 <- fh %>%
  ggplot(aes(x=V3)) +
  geom_histogram(aes(y=after_stat(cumsum(..count..))), binwidth=50)+
  stat_bin(aes(y=after_stat(cumsum(..count..))), bins=50,geom="line",color="blue") + 
  ylab("Number of bases") + 
  xlab("Coverage") + 
  theme_bw() + 
  ggtitle(args[1])


write.csv(file=args[3],
          threshold(fh),
          row.names=F)

ggsave(file=args[2],
       plot=plot_grid(p1,p2, ncol=2),
       dpi=600,
       width=7,
       height=5)
