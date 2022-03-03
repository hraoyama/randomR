# install.packages(c("igraph","statnet","network","igraphdata","intergraph", "RSiena", "RSienaTest"))
# install.packages(c("sna","circlize","visNetwork","networkD3"))
# install.packages(c("space","visNetwork","threejs","htmlwidgets","fields","LambertW"))

require(devtools)
install_github("DougLuke/UserNetR")
install_github("gastonstat/arcdiagram")
require(rgl)
require(arcdiagram)
require(UserNetR)

curve(1/sqrt(3*x) - 0.5*x, from=0, to=2, xlab="x", ylab="y"); grid()

library(igraph)
# g <- graph_from_literal(A--B, B-+C, C-+A)

g <- graph_from_literal(A--B, A--C, A--D, B--C, C--D, D--E)
transitivity(g, type = "global") # global clustering coefficient
3*2 / ncol(combn(5,3))
sum(transitivity(g, type = "local")[1:4])/5
g = graph_from_literal(A--C,A--B, B--C, B--D, D--E,E--F, E--G, F--G)
g_el = as_edgelist(g)
g_am = as_adjacency_matrix(g)
plot(g)
g2 = make_star(n=10,mode = "undirected")
plot(g2)
betweenness(g2)
betweenness(g)
library(igraphdata)
data("USairports")
graph_attr(USairports)
vertex_attr_names(USairports)
vertex_attr(USairports, "City")
plot(USairports)
