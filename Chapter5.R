### Chi - Square Contigency Analysis
rm(list=ls())
library(dplyr)
library(ggplot2)


setwd("C:/Users/14842/Documents/GettingStartedwR/datasets-master/datasets-master")
getwd()
ladybird<-read.csv("ladybirds_morph_colour.csv")
glimpse(ladybird)

total<- ladybird%>%
  group_by(Habitat,morph_colour)%>%
  summarise(total.number=sum(number))


### Summarize Raw Data with a Bar Chart

ggplot(total,aes(x=Habitat,y=total.number,fill=morph_colour))+
  geom_bar(stat='identity',position='dodge')

### Manually choose fill colors

ggplot(total,aes(x=Habitat,y=total.number,fill=morph_colour))+
  geom_bar(stat='identity',position='dodge')+scale_fill_manual(values=c(black="black",red="red"))


### Chi Square
## First we need to transform dataframe into a matrix

lady.mat<-xtabs(number~Habitat+morph_colour,data=ladybird)

chisq.test(lady.mat)
## This tells us to reject the NULL hypothesis
  #How do I know whether this is testing hypothesis or null?
lady.chi<-chisq.test(lady.mat)
names(lady.chi)
lady.chi$expected

#### Textbook

sample(1:40,5)
sample(1:40,5)
sample(40,5)
sample(8,5,replace=TRUE)
sample(c("H","T"),10,replace=T)
sample(c("success","fail"),10,replace=T,prob=c(0.9,0.1))
1/prod(40:36)
prod(5:1)/prod(40:36)
## Following compute 40!/5!35!
choose(40,5)
### probabiliy
1/choose(40,5)

## density function, used to draw bell curve of normal distribution

x<-seq(-4,4,0.1)
plot(x,dnorm(x))
plot(x,dnorm(x),type="l")

curve(dnorm(x),from=-4,to=4)

## for a discrete distribution, where variables can only take on distinct values, draw a pin diagram

x<-0:50
plot(x,dbinom(x,size=50,prob=.33),type="h")
## h stands for histogram