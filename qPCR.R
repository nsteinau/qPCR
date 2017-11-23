# Load packages we will need

library(reshape2)
library(readxl)
library(dplyr)
library(ggplot2)

# Read in the CT data from an xl file. Format needs to be of form Sample - Target - Value
# with only a one-line header and no extraneous text
f <- readline(prompt = "Enter the file name: ")
excel <- as.data.frame(read_excel(f))

# Reformat dataframe
dimnames(excel)[[2]] = c("sample","target","CT")
excel <- arrange(excel,target,sample,CT)
excel_gapdh<-filter(excel,target=='GAPDH')
excel <- filter(excel,target != 'GAPDH')
r<-nrow(excel)/nrow(excel_gapdh)
excel$GAPDH_CT <- rep(excel_gapdh$CT,r)
excel[excel == 'Undetermined'] <- 40
excel$target <- factor(excel$target)
excel$sample <- factor(excel$sample)
excel$CT <- as.numeric(excel$CT)
excel$GAPDH_CT <- as.numeric(excel$GAPDH_CT)

# Add new values using mutate
excel <- mutate(excel,CT_diff = CT - GAPDH_CT)
excel <- mutate(excel,power = 2 ^ (0-CT_diff))
std<-function(x){sd(x)/sqrt(length(x))}
a <- funs(mean,std)
excel_sum <- excel %>% group_by(sample,target) %>% summarise_at('power',a)

# If the qPCR data is from only two samples, this extra step will calculate the 
# two-sample t test p-value between the sample groups for each gene (target)
if(length(unique(excel$sample)) <= 2){
  pvals <- vector("list",r)
  pos <- 1
  for(i in unique(excel$target)){
    temp <- excel %>% filter(target == i)
    pvals[[pos]] <- t.test(temp$power[temp$sample == unique(excel$sample)[[1]]],temp$power[temp$sample == unique(excel$sample)[[2]]])
    pos <- pos + 1
  }
  remove(temp)
}

## This calculates error bars
limits <- aes(ymax = mean + std, ymin = mean - std)

## Set plot theme and plot
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = .5, face = "bold"), 
             axis.title = element_text(face = "bold"),
             axis.text = element_text(size = 10, face = "bold"),
             rect = element_rect(colour = "black", size = .75))

g <- ggplot(excel_sum, aes(x=target, y=mean, fill=sample)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(limits,
                width=.2,                    
                position=position_dodge(.9)) +
  labs(title = "Gene Expression", x = 'Target', y = 'Expression to GAPDH')

genes <- unique(excel$target)
plots <- vector("list", length = length(genes))
theme_update(aspect.ratio = 2)
pos <- 1
for(i in genes){
  plots[[pos]] <- ggplot(excel_sum[excel_sum$target == i,], aes(x=sample, y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity", fill = c("red","blue"),
             color = "black", size = .75) +
    geom_errorbar(limits,
                  width=.2,                    
                  position=position_dodge(.9), size = .75) +
    labs(title = paste(i,"( p =",signif(pvals[[pos]]$p.value, 2),")", sep = " "),x = 'Condition', y = 'Delta Ct / GAPDH')
  pos <- pos + 1
}

pos <- 1

## Save all plots to current working directory
for(i in plots){
  ggsave(filename = paste(as.character(genes[[pos]]),".png", sep = ""), plot = i, width = 3, height = 5.5)
  pos <- pos + 1
}

ggsave(filename = "combined.png", g, scale = 2)


