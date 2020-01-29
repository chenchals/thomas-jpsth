# https://www.r-graph-gallery.com/309-intro-to-hierarchical-edge-bundling.html
#
# Let’s start by creating the hierarchic structure with R. A hierarchic
# structure is basically a set of nodes, with edges linking nodes. We often
# accompany it with a second data frame that gives features for each node of the
# first data frame.
#
# Let’s build these 2 tables:

# createDataInput ---------------------------------------------------------
# Libraries
library(ggraph)
library(igraph)

# create a data frame giving the hierarchical structure of your individuals. 
# Origin on top, then groups, then subgroups
d1 <- data.frame(from="origin", to=paste("group", seq(1,10), sep=""))
# > dim(d1)
# [1] 10  2
# > head(d1)
# from     to
# 1 origin group1
# 2 origin group2
# 3 origin group3
# 4 origin group4
# 5 origin group5
# 6 origin group6
# > 

d2 <- data.frame(from=rep(d1$to, each=10), to=paste("subgroup", seq(1,100), sep="_"))
# > dim(d2)
# [1] 100   2
# > head(d2)
# from         to
# 1 group1 subgroup_1
# 2 group1 subgroup_2
# 3 group1 subgroup_3
# 4 group1 subgroup_4
# 5 group1 subgroup_5
# 6 group1 subgroup_6
# > 

hierarchy <- rbind(d1, d2)
# > dim(hierarchy)
# [1] 110   2
# > hierarchy[seq(1,40,by=5),]
# from          to
# 1  origin      group1
# 6  origin      group6
# 11 group1  subgroup_1
# 16 group1  subgroup_6
# 21 group2 subgroup_11
# 26 group2 subgroup_16
# 31 group3 subgroup_21
# 36 group3 subgroup_26
# > 

# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) ) 
# > dim(vertices)
# [1] 111   1
# > vertices[seq(1,30,by=5),]
# [1] origin      group5      group10     subgroup_5  subgroup_10 subgroup_15
# 111 Levels: group1 group10 group2 group3 group4 group5 group6 group7 group8 group9 origin subgroup_1 ... subgroup_99
# > 
#
# visualize ---------------------------------------------------------------
# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!
# > mygraph
# IGRAPH 9d10030 DN-- 111 110 -- 
#   + attr: name (v/c)
# + edges from 9d10030 (vertex names):
#   [1] origin->group1      origin->group2      origin->group3      origin->group4      origin->group5      origin->group6      origin->group7     
# [8] origin->group8      origin->group9      origin->group10     group1->subgroup_1  group1->subgroup_2  group1->subgroup_3  group1->subgroup_4 
# [15] group1->subgroup_5  group1->subgroup_6  group1->subgroup_7  group1->subgroup_8  group1->subgroup_9  group1->subgroup_10 group2->subgroup_11
# [22] group2->subgroup_12 group2->subgroup_13 group2->subgroup_14 group2->subgroup_15 group2->subgroup_16 group2->subgroup_17 group2->subgroup_18
# [29] group2->subgroup_19 group2->subgroup_20 group3->subgroup_21 group3->subgroup_22 group3->subgroup_23 group3->subgroup_24 group3->subgroup_25
# [36] group3->subgroup_26 group3->subgroup_27 group3->subgroup_28 group3->subgroup_29 group3->subgroup_30 group4->subgroup_31 group4->subgroup_32
# [43] group4->subgroup_33 group4->subgroup_34 group4->subgroup_35 group4->subgroup_36 group4->subgroup_37 group4->subgroup_38 group4->subgroup_39
# [50] group4->subgroup_40 group5->subgroup_41 group5->subgroup_42 group5->subgroup_43 group5->subgroup_44 group5->subgroup_45 group5->subgroup_46
# + ... omitted several edges
# > 
# With igraph: 
plot(mygraph, vertex.label="", edge.arrow.size=0, vertex.size=2)

# With ggraph:
ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_link() +
  theme_void()

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal() +
  theme_void()

# Add_a_few_connections -----------------------------------------------------
# Now, let’s add a second input to our data: connections. Suppose that nodes 18,
# 20 and 30 are connected to nodes 19, 50 and 70 respectively.
#
# An obvious solution to represent this link could be to add a straight line
# (left). The hierarchical edge bundling method does almost that. But it curves
# the lines to make thelm follow the edges of our structure (right).
#
# This method offers a tension parameters which controls how much we want to
# curve the lines.

# left: What happens if connections are represented with straight lines
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(alpha=0.1) +
  geom_conn_bundle(data = get_con(from = c(18,20,30), to = c(19, 50, 70)), alpha=1, width=1, colour="skyblue", tension = 0) +
  theme_void()

# right: using the bundle method (tension = 1)
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(alpha=0.1) +
  geom_conn_bundle(data = get_con(from = c(18,20,30), to = c(19, 50, 70)), alpha=1, width=1, colour="skyblue", tension = 1) +
  theme_void()

# Hierarchical_Edge_Bundling ------------------------------------------------
# Usually connections are stored in another data frame, here called connect. We
# have to pass it to ggraph to automatically plot all the connections. You get a
# hierarchical edge bundling chart.
#
# Note: ggraph expect nodes to be called following their id. Thus, it is
# necessary to get them using the match() function.

# create a dataframe with connection between leaves (individuals)
all_leaves <- paste("subgroup", seq(1,100), sep="_")
connect <- rbind( 
  data.frame( from=sample(all_leaves, 100, replace=T) , to=sample(all_leaves, 100, replace=T)), 
  data.frame( from=sample(head(all_leaves), 30, replace=T) , to=sample( tail(all_leaves), 30, replace=T)), 
  data.frame( from=sample(all_leaves[25:30], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)), 
  data.frame( from=sample(all_leaves[75:80], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)) 
)
# > dim(connect)
# [1] 190   2
# > head(connect)
# from          to
# 1 subgroup_62 subgroup_32
# 2  subgroup_8 subgroup_15
# 3 subgroup_44 subgroup_45
# 4 subgroup_36 subgroup_67
# 5 subgroup_82 subgroup_15
# 6  subgroup_3 subgroup_12
# > 

# The connection object must refer to the ids of the leaves:
from <- match( connect$from, vertices$name)
to <- match( connect$to, vertices$name)
# > head(from)
# [1] 73 19 55 47 93 14
# > head(to)
# [1] 43 26 56 78 26 23
# > 
# plot
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = 0) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()

# plot
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = 0.9) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()
# Conclusion ----------------------------------------------------------------
# This blogpost defined what hierarchical edge bundling is, and demonstrates how
# to build a basic one with R and ggraph. Now, go to the next level and learn
# how to customize this figure.




