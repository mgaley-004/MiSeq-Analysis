library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

print("This is an example of an R Scipt")
setwd("~")
mdf <- read.csv("example_metadata.csv", header=TRUE, row.names=1)
mdf.summary <- mdf %>% group_by(Species, Age, Sex, Reprod) %>% tally()
mdf.summary <- data.frame(mdf.summary)
print("I'm running a summary")
write.csv(mdf.summary, "mdfsummary.csv")

print("I'm printing a graphic")
theme_set(theme_classic(base_size=12, base_family="serif"))
ggplot(mdf.summary, aes(x=Species, y=n, fill=Sex)) + geom_bar(stat="identity") + scale_fill_manual(values=c("black", "gray"))

ggsave("mdftallies.png", width=5, height=5, units="in", dpi=300)

print("Script completed succesfully!")