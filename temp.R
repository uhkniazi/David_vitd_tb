# temp file 
# first look at data - 6/02/2015

source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGraph.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGeneAnnotation.R')

#library(Biobase)
library(DESeq2)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

mDat = as.matrix(read.csv(file.choose(), sep = ' ', header=T))
f = colnames(mDat)
f = gsub('^H._(\\w+)', '\\1', f)
f = factor(f)
dfDesign = data.frame(condition=f, row.names = colnames(mDat))

# call deseq2 constructor
oDseq = DESeqDataSetFromMatrix(mDat, dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)

# get the results for each comparison
oRes.rv.vs.med = results(oDseq, contrast = c('condition', 'Rv', 'med'))
oRes.rvd.vs.med = results(oDseq, contrast = c('condition', 'RvD', 'med'))
oRes.rv.vs.rvd = results(oDseq, contrast = c('condition', 'Rv', 'RvD'))

# get results with significant p-values
dfRv.vs.med = as.data.frame(oRes.rv.vs.med[which(oRes.rv.vs.med$padj < 0.1),])
dfRvD.vs.med = as.data.frame(oRes.rvd.vs.med[which(oRes.rvd.vs.med$padj < 0.1),])
dfRv.vs.RvD = as.data.frame(oRes.rv.vs.rvd[which(oRes.rv.vs.rvd$padj < 0.1),])

# create a incidence matrix with rows as genes and columns as types of comparisons
a = rownames(dfRv.vs.med)
b = rownames(dfRvD.vs.med)
c = rownames(dfRv.vs.RvD)
m = c(a, b, c)
m = unique(m)
mNames = matrix(F, nrow = length(m), ncol = 3, dimnames = list(m, c('a', 'b', 'c')))
f = m %in% a
mNames[,'a'] = f
f = m %in% b
mNames[,'b'] = f
f = m %in% c
mNames[,'c'] = f

# get the groups in the data
set.seed(1)
d = dist(mNames, method = 'binary')
hc = hclust(d)
t = cutree(hc, k = 6)
dfGroups = data.frame(mNames, group=t)
# get the ensembl names
cvEns = rownames(dfGroups)
# get the enterez names for these genes
dfEnt = f_dfEnsemblToEntrezID.org.Hs(cvID = cvEns)
# remove any duplicated data as versions of the gene are old in one annotation vs the newer enterez one
dfEnt = na.omit(dfEnt)
f = !duplicated(dfEnt$ENTREZID)
dfEnt = dfEnt[f,]
f = !duplicated(dfEnt$ENSEMBL)
dfEnt = dfEnt[f,]
# subset only these genes
dfGroups.full = dfGroups
f = rownames(dfGroups) %in% dfEnt$ENSEMBL
dfGroups = dfGroups[f,]
table(dfGroups$group)
# change rownames to enterez
m1 = rownames(dfGroups)
m2 = dfEnt$ENSEMBL
i = match(m1, m2)
# i is in mapping of m1 to m2
n = dfEnt[i, 'ENTREZID']
rownames(dfGroups) = n
dfGroups.full = dfGroups

# select the comparison of choice
dfGroups = dfGroups[dfGroups$group %in% c(1, 3, 5, 6),]
# get the log fold change for the chosen condition
m1 = rownames(dfGroups)
m2 = dfEnt$ENTREZID
i = match(m1, m2)
n = dfEnt[i, 'ENSEMBL']
dfGroups = cbind(dfGroups, dfEnt[i,])
temp = rbind(dfRvD.vs.med, dfRv.vs.med, dfRv.vs.RvD)
dfGroups = cbind(dfGroups, temp[dfGroups$ENSEMBL, ])
# convert the data to a graph format
library(igraph)
# get the dataframe to make bipartite graph using the CGeneAnnotation class
oGenes = CGeneAnnotation(rownames(dfGroups), org.Hs.eg.db)
# get the gene ontology to gene name mapping table
#dfGraph = CGeneAnnotation.dfgetGObyOntology(oGenes, 'BP')
dfGraph = CGeneAnnotation.dfgetReactome(oGenes)
dfGraph = na.omit(dfGraph)
# create bipartite graph
oIGbp = graph.data.frame(dfGraph, directed = F)
# set the type variable to make graph bipartite
#f = (grepl(pattern='GO:', V(oIGbp)$name))
f = rep(c(T, F), times = c(length(unique(dfGraph$ENTREZID)),length(unique(dfGraph$REACTOMEID))))
# V(oIGbp)$type = TRUE
# V(oIGbp)[f]$type = FALSE
# sanity check - is graph bipartite
V(oIGbp)$type = f
is.bipartite(oIGbp)

# create labels for the gene names, these short labels are used for plotting
f = V(oIGbp)$type
c = length(V(oIGbp)[f])
V(oIGbp)[f]$label = as.character(1:c)
# labels for go terms or reactome ids
c = length(V(oIGbp)[!f])
V(oIGbp)[!f]$label = paste('R', 1:c, sep='')
# assign a color, and shape to vertices
# lightblud for genes
V(oIGbp)[f]$color = 'lightblue'
# grey for go terms
V(oIGbp)[!f]$color = 'grey'
# assign a shape for each vertex
# circle for genes
V(oIGbp)[f]$shape = 'circle'
# rectangle for goterms
V(oIGbp)[!f]$shape = 'rectangle'
# V(oIGGoTableC)[!f]$size = 20
# V(oIGGoTableC)[f]$size = 30
# remove higher level go terms with many connections
ivDegGo = degree(oIGbp, V(oIGbp)[!f])
t = log(ivDegGo)
summary(t)
# which distribution can approximate the frequency of go terms
hist(t, prob=T, main='distribution of go terms with number associated genes', breaks=seq(-0.5, 8.5, by=1))
dn = dnbinom(0:8, size = mean(t), mu = mean(t))
dp = dpois(0:8, mean(t))
lines(0:8, dn, col='black', type='b')
lines(0:8, dp, col='red', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# try a cutoff of exp(1.5) based on the nbinom curve fit
# a poisson distribution with mean(t) fits will use this as cutoff
i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
c = names(which(ivDegGo>i))
v = V(oIGbp)[c]
oIGbp = delete.vertices(oIGbp, v)
# delete any orphan genes i.e. with no goterms
d = degree(oIGbp)
oIGbp = delete.vertices(oIGbp, which(d == 0))
# set the flag again, True is for goterm vertices and false for gene vertices
#f = (grepl(pattern='GO:', V(oIGbp)$name))
f = !V(oIGbp)$type
# save the names of associate goterms and gene vertex
m = get.edgelist(oIGbp)
lGene.to.GO = split(m[,2], m[,1])
mGene.to.GO = m

# get gene names into a table
cvVerGenes = V(oIGbp)[!f]$name
dfVerGenes = dfGroups[rownames(dfGroups) %in% cvVerGenes,]
lab = V(oIGbp)[rownames(dfVerGenes)]$label
dfVerGenes$label = lab
# assign group labels to vertex properties
V(oIGbp)[rownames(dfVerGenes)]$group = dfVerGenes$group
V(oIGbp)[rownames(dfVerGenes)]$log2FoldChange = dfVerGenes$log2FoldChange

# assign symbols to genes
n = V(oIGbp)$name
dfSym = select(org.Hs.eg.db, n, 'SYMBOL', 'ENTREZID')
V(oIGbp)$symbol = dfSym$SYMBOL

########### analysis
######### graph projection
# create the CGraph object and calculate obs to exp weights after projection
obj = CGraph(oIGbp)

# create a projection of the graph 
oIGProj = CGraph.getProjectedGraph(obj)#bipartite.projection(oIGbp, which='TRUE')
# remove those genes with no connections i.e. degrees == 0
ivDegGo = degree(oIGProj)
# get the names of those vertices 
c = names(which(ivDegGo==0))
v = V(oIGProj)[c]
oIGProj = delete.vertices(oIGProj, v)
# show overexpressed genes as circle and underexpressed as square
V(oIGProj)$shape = ifelse(V(oIGProj)$log2FoldChange > 0, 'circle', 'square')
# switch the weights with obs to exp ratio
E(oIGProj)$weight_old = E(oIGProj)$weight
E(oIGProj)$weight = E(oIGProj)$ob_to_ex
x = E(oIGProj)$ob_to_ex
# set minimum value to 1 i.e. below 1 to 1
# x[x < 1] = 1
# E(oIGProj)$weight = 10*log(x)
# delete edges with higher than 0.05 quantile under normal distribution
x = log(x)
hist(x, prob=T, breaks = seq(min(x)-0.5, max(x)+0.5, 1))
dn = dnorm(0:9, mean(x), sd(x))
lines(0:9, dn, col='red', type='b')
lines(density(x), col='blue')
i = qnorm(0.95, mean(x), sd(x), lower.tail = T)
# or use a quantile 
i = quantile(x, 0.75)
i = which(x <= i)
# remove these edges
oIGProj = delete.edges(oIGProj, edges = i)
#E(oIGProj)$weight = E(oIGProj)$ob_to_ex
# delete the orphan vertices again
ivDegGo = degree(oIGProj)
# get the names of those vertices 
c = names(which(ivDegGo==0))
v = V(oIGProj)[c]
oIGProj = delete.vertices(oIGProj, v)
p.old = par(mar=c(1,1,1,1))

plot(oIGProj, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(oIGProj)$group, vertex.frame.color=NA)
legend('topright', legend = unique(V(oIGProj)$group), fill=unique(V(oIGProj)$group))

# get the count data from the deseq object
# and create a correlation matrix
#mCounts = assay(rlog(oDseq))
mCounts = counts(oDseq, normalized=T)
mCounts = na.omit(log(mCounts))
f = is.finite(rowSums(mCounts))
mCounts = mCounts[f,]
# data quality not right
plot(density(rowMeans(mCounts)))

# get the names of the genes in enterez ids from the graph
n = V(oIGProj)$name
# subset the mapping dataframe on those ids to get matching ensembl ids
n2 = dfEnt[dfEnt$ENTREZID %in% n, 'ENSEMBL' ]
# subset the count matrix choosing only these ensembl ids
mCounts.sub = mCounts[rownames(mCounts) %in% n2,]
# replace these ensembl ids with matching i.e. correspoinding enterez ids
n = rownames(mCounts.sub)
# find these matching rownames in the mapping dataframe
n2 = dfEnt[dfEnt$ENSEMBL %in% n, ]
i = match(n, n2$ENSEMBL)
# i preserves the order of match - replace old names by these new names
n = n2$ENTREZID[i]
rownames(mCounts.sub) = n
# create the incidecne matrix
par(p.old)
d = dist(mCounts.sub)
hc = hclust(d)
plot(hc)
m = cor(t(mCounts.sub))
par(p.old)
diag(m) = 0
hist(m)
# create the graph
oIGcor = graph.adjacency(m, mode='min', weighted=T)
# check the distribution of the correlations
c = E(oIGcor)$weight
hist(c)
summary(c)
E(oIGcor)$cor = E(oIGcor)$weight
t = abs(c[c < 0])
hist(t, prob=T)
dn = dnorm(seq(0, 1, 0.01), mean(t), sd(t))
lines(seq(0, 1, 0.01), dn)
summary(t)
dn = dnorm(seq(0, 1, 0.01), median(t), mad(t))
lines(seq(0, 1, 0.01), dn, col='red')
# choose around 0.8
f = which(c <=0.8) #f = which((c >= -0.7 & c <= 0.7))
oIGcor = delete.edges(oIGcor, edges = f)

# intersect the 2 graphs
l = list(oIGProj, oIGcor)
ig.1 = graph.intersection(l)
E(ig.1)$weight = E(ig.1)$weight_1
par(mar=c(1,1,1,1))
# color the positive and negative edges differently
col = ifelse(E(ig.1)$weight_2 < 0, 'red', 'black')
E(ig.1)$color = col
plot(ig.1, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(ig.1)$group, edge.width=0.5,
     vertex.frame.color=NA)
legend('topright', legend = unique(V(ig.1)$group), fill=unique(V(ig.1)$group))
# look for clusters i.e. connected components in the graph
# remove genes with no connections
d = degree(ig.1)
#par(p.old)
# hist(d, prob=T, breaks=seq(-0.5, 35.5, by=1))
# summary(d)
# dn = dnbinom(0:12, size=mean(d), mu=mean(d))
# dp = dpois(0:12, mean(d))
# legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# lines(0:12, dn, col='black', type='b')
# lines(0:12, dp, col='red', type='b')
# # choose a cutoff
# qpois(0.05, lambda = mean(d), lower.tail = F)
# qnbinom(0.05, size = mean(d), mu = mean(d), lower.tail = F)
d = which(degree(ig.1) == 0)
ig.2 = delete.vertices(ig.1, V(ig.1)[d])
par(mar=c(1,1,1,1))
plot(ig.2, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(ig.2)$group, edge.width=0.5)
legend('topright', legend = 1:6, fill=1:6)
# which are the largest components
cl = clusters(ig.2)
#table(cl$csize)
i = which(cl$csize == 9)
# get the index of these vertex ids # choose one cluster at a time for community checking
f = which(cl$membership %in% i)
# OR use this
f = which(cl$membership == 1)
ig.3 = induced.subgraph(ig.2, vids = f)
par(mar=c(1,1,1,1))
d = degree(ig.3)
d = log(d)
plot(ig.3, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(ig.3)$group, edge.width=0.5,
     vertex.frame.color=NA)
legend('topright', legend = unique(V(ig.3)$group), fill=unique(V(ig.3)$group))
# find the communities within this graph
# assign a new weight for community as zero is not acceptable
# E(ig.3)$weight = E(ig.3)$weight+1
com = edge.betweenness.community(ig.3)
com = fastgreedy.community(ig.3)
plot(com, ig.3, vertex.label=V(ig.3)$label, 
     vertex.size=1, 
     layout=layout.fruchterman.reingold, 
     edge.width=0.5)
write.graph(ig.3, 'Data_export/group_b_largest_cluster.graphml', format = 'graphml')
##### community analysis
ivCom = membership(com)
dfCom = data.frame(id=V(ig.3)$name, symbol=V(ig.3)$symbol, lfc=V(ig.3)$log2FoldChange, com=ivCom)
par(p.old)
hist(dfCom$lfc)
t.test(dfCom$lfc)
dfCom$lfc.c = scale(dfCom$lfc)
t.test(dfCom$lfc.c)
f = sizes(com)
f = f[f > 2]
i = names(f)
i = as.numeric(i)
dfCom.sub = dfCom[dfCom$com %in% i,]
dim(dfCom.sub)
com.pval = tapply(dfCom.sub$lfc, INDEX = dfCom.sub$com, FUN = function(x) t.test(x)$p.value)
com.sig = names(which(com.pval < 0.05))
com.sig = as.numeric(com.sig)
temp = dfCom.sub[dfCom.sub$com %in% com.sig,]
temp = temp[order(temp$com),]

# 
####### betweenness analysis

ivBet = betweenness(ig.3, directed = F)
ivBet['240']
summary(ivBet)
temp = ivBet
temp = temp[temp >= quantile(temp, 0.50)]
n = names(temp)
length(n)
head(n)
ig.4 = induced.subgraph(ig.3, vids = n)

par(mar=c(1,1,1,1))
plot(ig.4, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(ig.4)$group, edge.width=0.5,
     vertex.frame.color=NA)
legend('topright', legend = unique(V(ig.4)$group), fill=unique(V(ig.4)$group))
write.graph(ig.4, file = 'Data_export/group_b_betweenness.graphml', format = 'graphml')




c = rainbow(length(com))
w = E(ig.3)$weight
par(oma=c(0,0,0,0), mar=c(0.2, 0.2, 0.2, 0.2))
plot(ig.3, vertex.size=d, 
     vertex.label.cex=d/max(d)+0.1, 
     #vertex.label=V(ig.3)$name,
     layout=layout.fruchterman.reingold,#(ig.3, weight=(E(ig.3)$weight)), 
      edge.width=1, edge.curved=F, #vertex.frame.color=NA,#vertex.frame.color=c[membership(com)],
     vertex.color=V(ig.3)$group, vertex.label.dist=0)#,
     #edge.label=round(E(ig.3)$weight,0), edge.label.color='red', edge.label.cex=w/max(w)+0.1)
legend('topright', legend = unique(V(ig.3)$group), fill=unique(V(ig.3)$group))
# for tkplot
ig.4 = ig.3
V(ig.4)$color = ifelse(V(ig.4)$log2FoldChange > 0, 'lightblue', 'red')
#V(ig.4)['240']$color = 'green'
tkplot(ig.4, #vertex.size=d, 
     #vertex.label.cex=d/max(d)+0.1, 
     vertex.label=V(ig.4)$name,
     layout=layout.fruchterman.reingold,#(ig.3, weight=(E(ig.3)$weight)), 
     edge.curved=T)#, vertex.color=V(ig.4)$group + 1)#, #vertex.frame.color=NA,#vertex.frame.color=c[membership(com)],
sizes(com)
length(membership(com))
# genes in each community
# which go terms are more common in each community
s = unique(membership(com))
dfGeneLabels = NULL
dfGO = NULL
for (i in 1:length(s)){
  # get the index number of members of current community i.e. i
  c = which(membership(com) == s[i])
  # get the enterez id names for these genes
  n = V(ig.3)[c]$name
  # get the label numbers for these genes
  lab = dfVerGenes[rownames(dfVerGenes) %in% n,'label']
  # get the gene annotation for these genes
  ob = CGeneAnnotation(n, org.Hs.eg.db)
  # get the gene names from the CGeneAnnotation object and create dataframe
  df = data.frame(CGeneAnnotation.dfGetGeneName(ob), n, lab, community=s[i])
  # add information to dataframe
  dfGeneLabels = rbind(dfGeneLabels, df) 
  # select the go ids for these genes
  df = sapply(seq_along(n), function(x) lGene.to.GO[[n[x]]])
  df = unlist(df)
  #select(org.Hs.eg.db, keys=n, columns='GO', keytype='ENTREZID')
  #df = df[df$ONTOLOGY == 'BP',]
  # sort the go ids based on their frequencies
  temp = as.data.frame(table(df))
  temp = temp[order(temp$Freq, decreasing = T),]
  # select the top 2 go ids
  go = temp$df[1:5]
  # get the description of these go ids
  #temp = select(GO.db, keys = as.character(go), columns = c('ONTOLOGY', 'DEFINITION'), keytype = 'GOID')
  temp = select(reactome.db, keys=as.character(go), columns=c("PATHNAME", "REACTOMEID"), keytype = 'REACTOMEID')
  temp$community = s[i]
  dfGO = rbind(dfGO, temp)
}
# i = 37
# i = which(membership(com) == i)
# n = V(ig.3)[i]$name
# select(org.Hs.eg.db, keys=n, columns='GENENAME', keytype='ENTREZID')
# df = select(org.Hs.eg.db, keys=n, columns='GO', keytype='ENTREZID')
# 
# temp = induced.subgraph(ig.3, i)
# plot(temp, vertex.label=V(temp)$name)
# plot(temp, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.color=V(temp)$group, edge.width=0.5)





cvGenes.C = rownames(dfRv.vs.RvD)
cvGenes.C = f_dfEnsemblToEntrezID.org.Hs(cvID = cvGenes.C)
cvGenes.C = cvGenes.C$ENTREZID
cvGenes.C = cvGenes.C[!is.na(cvGenes.C)]
lGenes.C = lapply(seq_along(lTemp), function(x) CGeneAnnotation(cvGenes.C[x], org.Hs.eg.db))
names(lGenes.C) = cvGenes.C
dfGO = NULL

for (i in 1:length(lGenes.C)){
  dfGO = rbind(dfGO, CGeneAnnotation.dfgetGObyOntology(lGenes.C[[i]]))
}

cvGO.C = unique(dfGO$GO)
# get the goterms for other 3 classes


# #### testing graph weights methods
# f = read.csv(file.choose(), header = T, sep = ',', row.names=1)
# g = graph.incidence(f)
# vertex.attributes(g)
# V(g)$type = !(V(g)$type)
# is.bipartite(g)
# vertex.attributes(g)
# plot(g, layout=layout.bipartite)
# obj = CGraph(g)
# plot(obj@ig, layout=layout.bipartite)
# 
# # assign probabilities to vertex of first kind
# CGraph.assign.marginal.probabilities = function(obj){
#   # vertex of the first kind will be assigned probabilities
#   # based on their relations with the vertices of the second kind
#   # flag to identify vertex types
#   f = V(obj@ig)$type
#   d = degree(obj@ig)
#   d = d[f]
#   # r is the total numbers of vertices of the second kind
#   r = sum(!f)
#   p = d/r
#   V(obj@ig)[f]$prob_marginal = p
#   obj@r = r
#   obj@f = f
#   return(obj)
# }
# 
# 
# CGraph.project = function(obj){
#   # project the graph in one dimension and
#   # assign weights based on observed to expected ratios
#   g.p = bipartite.projection(obj@ig, which = 'TRUE')
#   # get the matrix with rows representing each edge
#   m = get.edgelist(g.p)
#   w = E(g.p)$weight
#   # calculate observed ratio
#   # weight / r
#   ob = w / obj@r
#   # calculate expected 
#   mExp = cbind(V(g.p)[m[,1]]$prob_marginal, V(g.p)[m[,2]]$prob_marginal)
#   ex = mExp[,1] * mExp[,2]
#   E(g.p)$observed = ob
#   E(g.p)$expected = ex
#   E(g.p)$ob_to_ex = ob / ex
#   obj@ig.p = g.p
#   return(obj)
# }

library("RColorBrewer", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.1")
# generate some random data
d = runif(10, min = 1, max = 10)
# sort it in ascending order
d = sort(d)
# make a matrix with 2 rows , second row is same data in ascending order
mat = rbind(d, rev(d))
### generating your own colours
# get the colour code e.g. red
col2rgb('red')
# using the 3 values you get, generate the hex code for the colour
r = rgb(255, 0, 0, maxColorValue = 255)
y = rgb(255, 255, 0, maxColorValue = 255)
g = rgb(0, 255, 0, maxColorValue = 255)
# put colour in the palette
col = c(g, y, r)
# make heat map
heatmap(mat, col=col, Colv = NA)
# use reverse order of colours
heatmap(mat, col=rev(col), Colv=NA)

## use built in palette, using 4 colours
col = brewer.pal(4, 'RdYlGn')
heatmap(mat, col=col, Colv = NA)
# use reverse order of colours
heatmap(mat, col=rev(col), Colv=NA)

