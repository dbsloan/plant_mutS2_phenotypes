#Set working directory

#Designate which packages to load and library
#This is a long list of packages I use frequently but I don't think all of them are at use in this script.
#If any are causing problems, try removing them from the list.
package_list<-c("dplyr", "ape", "treeio", "ggtree", "tidytree", "seqinr", "ggplot2", "phytools", "phangorn", "viridis")

#Loop to check if package is installed and 'libraried'. If it isn't already, then install and library
for(k in 1:length(package_list)){
  library(package_list[k], character.only=TRUE)
}


#Read in the  tree (as a typical tree object)
tree<-ggtree::read.tree(file = "../RAxML_bipartitions.nwk")

ogroupNames <- c("Bacillus_subtilis_subsp_CAB14818.1", "Borreliella_burgdorferi_AAC66481.1",
            "Helicobacter_pylori_AAD07685.1", "Aquifex_aeolicus_AAC07247.1")

#root tree
tree_unrooted <- ape::unroot(tree)
tree_rooted <- ape::root(tree_unrooted, outgroup = ogroupNames, 
                         node = (mrca.phylo(tree_unrooted, ogroupNames)), 
                         resolve.root=TRUE)


#Read in domain data for the tree
#This file was generated using the web-based tool at: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
#then remove the first 7 lines 
domain_dat<-read.table(file = "domains.txt", sep = "\t", header = TRUE)

#Change first column name
names(domain_dat)[1]<-"Newick_label"

#Get sequence file
seqs<-seqinr::read.fasta(file = "seqs.fasta")

#add the sequence length data to the df
domain_dat_full<-right_join(domain_dat, data.frame(Newick_label=names(seqs), Seq_ln=getLength(seqs)))

#Change the classes in the data frame
domain_dat_full[,1]<-paste(domain_dat_full[,1]) #sequence names
domain_dat_full[,4]<-as.numeric(paste(domain_dat_full[,4])) #from
domain_dat_full[,5]<-as.numeric(paste(domain_dat_full[,5])) #to
domain_dat_full[,6]<-as.numeric(paste(domain_dat_full[,6])) #e-val
domain_dat_full[,7]<-as.numeric(paste(domain_dat_full[,7])) #bitscore
domain_dat_full[,8]<-paste(domain_dat_full[,8]) #accession ID
domain_dat_full[,9]<-paste(domain_dat_full[,9]) #short name
domain_dat_full[,10]<-paste(domain_dat_full[,10]) #incomplete
domain_dat_full[,11]<-paste(domain_dat_full[,11]) #superfamily
domain_dat_full[,12]<-as.numeric(paste(domain_dat_full[,12]))
#Make a new column that's the same as newick labels
domain_dat_full[,13]<-paste(domain_dat_full[,1])
names(domain_dat_full)[13]<-"TipLabels"

#Make a ggtree object - does NOT keep BS vals
p<-ggtree(tree_rooted, branch.length ='none', ladderize = TRUE) + scale_color_viridis(discrete = TRUE)

#Add tip names in as a facet
p2<-facet_plot(p, panel='Species Names',
               data=domain_dat_full, geom=geom_text, hjust = 0,
               mapping=aes(x=0, label= TipLabels), size=3)

#add seq length line
p4<-facet_plot(p2, panel = "domains", data = domain_dat_full, geom= geom_segment, 
               mapping = aes(x=0, xend=Seq_ln, y=y, yend=y), size=0.5, color='black')

#Add domains
p5<-facet_plot(p4, panel = "domains", data = domain_dat_full, geom=geom_segment, 
               aes(x=From, xend=To, y=y, yend=y, col=Short.name), size=3) +
  theme(legend.position = "right")

#Add names of domains (optional)
p6<-facet_plot(p5, panel='domains',
                data=domain_dat_full, geom=geom_text,
                mapping=aes(x=From+20, label= Short.name), size=2)

#Save tree as pdf
ggplot2::ggsave(filename = "Domains.pdf", 
                plot = p5,
                width = 8, 
                height = (length(tree$tip.label)/4)) 
