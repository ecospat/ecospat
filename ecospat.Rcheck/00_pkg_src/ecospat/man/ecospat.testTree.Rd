\name{ecospat.testTree}

\docType{data}

\alias{ecospat.testTree}

\title{Test Tree For The Ecospat package}

\description{The tree object is a phylogenetic tree of class 'phylo' (see read.tree) that contains data of 50 angiosperm species from the Western Swiss Alps.}

\format{ecospat.testTree is a tree contains the following species:

 [1] "Rumex_acetosa"            
 [2] "Polygonum_bistorta"       
 [3] "Polygonum_viviparum"      
 [4] "Rumex_crispus"            
 [5] "Cerastium_latifolium"     
 [6] "Silene_acaulis"           
 [7] "Gypsophila_repens"        
 [8] "Vaccinium_gaultherioides" 
 [9] "Soldanella_alpina"        
[10] "Cruciata_laevipes"        
[11] "Galium_album"             
[12] "Galium_anisophyllon"      
[13] "Galium_megalospermum"     
[14] "Gentiana_verna"           
[15] "Gentiana_bavarica"        
[16] "Gentiana_purpurea"        
[17] "Gentiana_lutea"           
[18] "Bartsia_alpina"           
[19] "Rhinanthus_alectorolophus"
[20] "Prunella_grandiflora"     
[21] "Acinos_alpinus"           
[22] "Plantago_alpina"          
[23] "Plantago_lanceolata"      
[24] "Veronica_officinalis"     
[25] "Veronica_aphylla"         
[26] "Veronica_alpina"          
[27] "Veronica_chamaedrys"      
[28] "Veronica_persica"         
[29] "Globularia_cordifolia"    
[30] "Globularia_nudicaulis"    
[31] "Myosotis_alpestris"       
[32] "Myosotis_arvensis"        
[33] "Aposeris_foetida"         
[34] "Centaurea_montana"        
[35] "Hieracium_lactucella"     
[36] "Leontodon_helveticus"     
[37] "Leontodon_autumnalis"     
[38] "Hypochaeris_radicata"     
[39] "Achillea_atrata"          
[40] "Achillea_millefolium"     
[41] "Homogyne_alpina"          
[42] "Senecio_doronicum"        
[43] "Adenostyles_glabra"       
[44] "Arnica_montana"           
[45] "Aster_bellidiastrum"      
[46] "Bellis_perennis"          
[47] "Doronicum_grandiflorum"   
[48] "Phyteuma_orbiculare"      
[49] "Phyteuma_spicatum"        
[50] "Campanula_rotundifolia" 

}


\references{Ndiribe, C., L. Pellissier, S. Antonelli, A. Dubuis, J. Pottier, P. Vittoz, A. Guisan and N. Salamin. 2013. Phylogenetic plant community structure along elevation is lineage specific. \emph{Ecology and Evolution}, \bold{3}, 4925-4939.
}

\author{Charlotte Ndiribe \email{charlotte.ndiribe@unil.ch}, Nicolas Salamin \email{nicolas.salamin@unil.ch} and Antoine Guisan \email{antoine.guisan@unil.ch}

}


\examples{
fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
tree <- read.tree(fpath)
plot(tree)
}