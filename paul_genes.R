#NAME: paul_genes.R
#DESC: input the gene lists from the files and create graphs
#AUTH: u.niazi@imperial.ac.uk
#Date: 20/03/2015


## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGraph.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGeneAnnotation.R')
library(igraph)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)



## data import and formatting
dfGrp.A = read.csv(file.choose(), header = T, sep='\t')
dfGrp.A$Group_no = as.factor(dfGrp.A$Group_no)

dfGrp.B = read.csv(file.choose(), header = T, sep='\t')
dfGrp.B$Group_no = as.factor(dfGrp.B$Group_no)
# remove duplicated genes
dfGrp.A = dfGrp.A[!duplicated(dfGrp.A$EntrezGene_ID),]
dfGrp.B = dfGrp.B[!duplicated(dfGrp.B$EntrezGene_ID),]

## data processing

## create graphs
# note: will give warnings for multiple key mappings
oGAgrp.A = CGeneAnnotation(as.character(dfGrp.A$EntrezGene_ID), db.obj = org.Hs.eg.db)
oGAgrp.B = CGeneAnnotation(as.character(dfGrp.B$EntrezGene_ID), db.obj = org.Hs.eg.db)
# create data frames of reactome ids mapped to gene ids to make bipartite graph
dfGraph.A = CGeneAnnotation.dfgetReactome(oGAgrp.A)
dfGraph.A = na.omit(dfGraph.A)
dfGraph.B = CGeneAnnotation.dfgetReactome(oGAgrp.B)
dfGraph.B = na.omit(dfGraph.B)
# create bipartite graph
oIGbp.A = graph.data.frame(dfGraph.A, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph.A$ENTREZID)),length(unique(dfGraph.A$REACTOMEID))))
V(oIGbp.A)$type = f
# sanity check - is graph bipartite
is.bipartite(oIGbp.A)
# for group B
oIGbp.B = graph.data.frame(dfGraph.B, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph.B$ENTREZID)),length(unique(dfGraph.B$REACTOMEID))))
V(oIGbp.B)$type = f
# sanity check - is graph bipartite
is.bipartite(oIGbp.B)

## graph cleaning
# remove high level reactome terms, as they create too many edges
f = V(oIGbp.A)$type
ivDegGo = degree(oIGbp.A, V(oIGbp.A)[!f])
t = log(ivDegGo)
summary(t)
# which distribution can approximate the frequency of reactome terms
hist(t, prob=T, main='distribution of reactome terms with number associated genes', breaks=seq(-0.5, 6.5, by=1))
# try negative binomial and poisson distributions
# parameterized on the means
dn = dnbinom(0:7, size = mean(t), mu = mean(t))
dp = dpois(0:7, mean(t))
lines(0:7, dn, col='black', type='b')
lines(0:7, dp, col='red', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# a poisson distribution with mean(t) fits well - use this as cutoff
i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
c = names(which(ivDegGo>i))
v = V(oIGbp.A)[c]
oIGbp.A = delete.vertices(oIGbp.A, v)

# repeat for graph of group B
f = V(oIGbp.B)$type
ivDegGo = degree(oIGbp.B, V(oIGbp.B)[!f])
t = log(ivDegGo)
summary(t)
# which distribution can approximate the frequency of reactome terms
hist(t, prob=T, main='distribution of reactome terms with number associated genes', breaks=seq(-0.5, 6.5, by=1))
# try negative binomial and poisson distributions
# parameterized on the means
dn = dnbinom(0:7, size = mean(t), mu = mean(t))
dp = dpois(0:7, mean(t))
lines(0:7, dn, col='black', type='b')
lines(0:7, dp, col='red', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# a poisson distribution with mean(t) fits well - use this as cutoff
i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
c = names(which(ivDegGo>i))
v = V(oIGbp.B)[c]
oIGbp.B = delete.vertices(oIGbp.B, v)

# delete any orphan genes left behind, i.e. genes with no reactome terms
d = degree(oIGbp.A)
oIGbp.A = delete.vertices(oIGbp.A, which(d == 0))
d = degree(oIGbp.B)
oIGbp.B = delete.vertices(oIGbp.B, which(d == 0))

# assign associated meta data to the gene vertices
rownames(dfGrp.A) = as.character(dfGrp.A$EntrezGene_ID)
rownames(dfGrp.B) = as.character(dfGrp.B$EntrezGene_ID)
f = V(oIGbp.A)$type
n = V(oIGbp.A)[f]$name
V(oIGbp.A)[n]$logFC = dfGrp.A[n,'logFC']
V(oIGbp.A)[n]$group = dfGrp.A[n,'Group_no']
dfSym = select(org.Hs.eg.db, n, 'SYMBOL', 'ENTREZID')
rownames(dfSym) = dfSym$ENTREZID
V(oIGbp.A)[n]$symbol = dfSym[n,'SYMBOL']
# repeat for group B
f = V(oIGbp.B)$type
n = V(oIGbp.B)[f]$name
V(oIGbp.B)[n]$logFC = dfGrp.B[n,'logFC']
V(oIGbp.B)[n]$group = dfGrp.B[n,'Group_no']
dfSym = select(org.Hs.eg.db, n, 'SYMBOL', 'ENTREZID')
rownames(dfSym) = dfSym$ENTREZID
V(oIGbp.B)[n]$symbol = dfSym[n,'SYMBOL']
rm(dfSym)

## graph projection to one dimension
# create the CGraph object and calculate obs to exp weights after projection
obj = CGraph(oIGbp.A)
# create a projection of the graph 
oIGProj.A = CGraph.getProjectedGraph(obj)
# repeat for group B
obj = CGraph(oIGbp.B)
# create a projection of the graph 
oIGProj.B = CGraph.getProjectedGraph(obj)
rm(obj)

## some genes are orphans as they don't share reactome terms with others, remove those
d = degree(oIGProj.A)
oIGProj.A = delete.vertices(oIGProj.A, which(d == 0))
d = degree(oIGProj.B)
oIGProj.B = delete.vertices(oIGProj.B, which(d == 0))

## house keeping
# show overexpressed genes as circle and underexpressed as square
V(oIGProj.A)$shape = ifelse(V(oIGProj.A)$logFC > 0, 'circle', 'square')
V(oIGProj.B)$shape = ifelse(V(oIGProj.B)$logFC > 0, 'circle', 'square')
# switch the weights with obs to exp ratio
E(oIGProj.A)$weight_old = E(oIGProj.A)$weight
E(oIGProj.A)$weight = E(oIGProj.A)$ob_to_ex
E(oIGProj.B)$weight_old = E(oIGProj.B)$weight
E(oIGProj.B)$weight = E(oIGProj.B)$ob_to_ex

## remove low observed to expected probabilities
x = E(oIGProj.A)$weight
# NOTE: this cutoff can be changed, the lower it is the more edges in the graph
i = quantile(x, 0.75)
i = which(x <= i)
# remove these edges
oIGProj.A = delete.edges(oIGProj.A, edges = i)
# repeat for group B
x = E(oIGProj.B)$weight
# NOTE: this cutoff can be changed, the lower it is the more edges in the graph
i = quantile(x, 0.75)
i = which(x <= i)
# remove these edges
oIGProj.B = delete.edges(oIGProj.B, edges = i)
# remove any orphan vertices 
d = degree(oIGProj.A)
oIGProj.A = delete.vertices(oIGProj.A, which(d == 0))
d = degree(oIGProj.B)
oIGProj.B = delete.vertices(oIGProj.B, which(d == 0))

## do a plot to visualize graph structure
p.old = par(mar=c(1,1,1,1))

plot(oIGProj.A, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, 
     vertex.color=V(oIGProj.A)$group, vertex.frame.color=NA)
legend('topright', legend = unique(V(oIGProj.A)$group), fill=unique(V(oIGProj.A)$group))
# graph for group B
plot(oIGProj.B, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, 
     vertex.color=V(oIGProj.B)$group, vertex.frame.color=NA)
legend('topright', legend = unique(V(oIGProj.B)$group), fill=unique(V(oIGProj.B)$group))

## export the graphs for use in cytoscape
write.graph(oIGProj.A, 'Data_export/groupA.graphml', format = 'graphml')
write.graph(oIGProj.B, 'Data_export/groupB.graphml', format = 'graphml')



