\name{ecospat.testData}

\docType{data}

\alias{ecospat.testData}

\title{Test Data For The Ecospat package}

\description{Data frame that contains vegetation plots data: presence records of 50 species, a set of environmental variables (topo-climatic) and SDM predictions for some species in the Western Swiss Alps (Canton de Vaud, Switzerland).
}

\usage{data("ecospat.testData")}

\format{A data frame with 300 observations on the following 96 variables.
  \describe{
    \item{\code{numplots}}{Number of the vegetation plot.}
    \item{\code{long}}{Longitude, in Swiss plane coordinate system of the vegetation plot.}
    \item{\code{lat}}{Latitude, in Swiss plane coordinate system of the vegetation plot.}
    \item{\code{ddeg}}{Growing degree days (with a 0 degrees Celsius threshold).}
    \item{\code{mind}}{Moisture index over the growing season (average values for June to August in mm day-1).}
    \item{\code{srad}}{The annual sum of radiation (in kJ m-2 year-1).}
    \item{\code{slp}}{Slope (in degrees) calculated from the DEM25.}
    \item{\code{topo}}{Topographic position (an integrated and unitless measure of topographic exposure.}
    \item{\code{Achillea_atrata}}{}
    \item{\code{Achillea_millefolium}}{}
    \item{\code{Acinos_alpinus}}{}
    \item{\code{Adenostyles_glabra}}{}
    \item{\code{Aposeris_foetida}}{}
    \item{\code{Arnica_montana}}{}
    \item{\code{Aster_bellidiastrum}}{}
    \item{\code{Bartsia_alpina}}{}
    \item{\code{Bellis_perennis}}{}
    \item{\code{Campanula_rotundifolia}}{}
    \item{\code{Centaurea_montana}}{}
    \item{\code{Cerastium_latifolium}}{}
    \item{\code{Cruciata_laevipes}}{}
    \item{\code{Doronicum_grandiflorum}}{}
    \item{\code{Galium_album}}{}
    \item{\code{Galium_anisophyllon}}{}
    \item{\code{Galium_megalospermum}}{}
    \item{\code{Gentiana_bavarica}}{}
    \item{\code{Gentiana_lutea}}{}
    \item{\code{Gentiana_purpurea}}{}
    \item{\code{Gentiana_verna}}{}
    \item{\code{Globularia_cordifolia}}{}
    \item{\code{Globularia_nudicaulis}}{}
    \item{\code{Gypsophila_repens}}{}
    \item{\code{Hieracium_lactucella}}{}
    \item{\code{Homogyne_alpina}}{}
    \item{\code{Hypochaeris_radicata}}{}
    \item{\code{Leontodon_autumnalis}}{}
    \item{\code{Leontodon_helveticus}}{}
    \item{\code{Myosotis_alpestris}}{}
    \item{\code{Myosotis_arvensis}}{}
    \item{\code{Phyteuma_orbiculare}}{}
    \item{\code{Phyteuma_spicatum}}{}
    \item{\code{Plantago_alpina}}{}
    \item{\code{Plantago_lanceolata}}{}
    \item{\code{Polygonum_bistorta}}{}
    \item{\code{Polygonum_viviparum}}{}
    \item{\code{Prunella_grandiflora}}{}
    \item{\code{Rhinanthus_alectorolophus}}{}
    \item{\code{Rumex_acetosa}}{}
    \item{\code{Rumex_crispus}}{} 
    \item{\code{Vaccinium_gaultherioides}}{}
    \item{\code{Veronica_alpina}}{}
    \item{\code{Veronica_aphylla}}{}
    \item{\code{Agrostis_capillaris}}{}
    \item{\code{Bromus_erectus_sstr}}{}
    \item{\code{Campanula_scheuchzeri}}{}
    \item{\code{Carex_sempervirens}}{}
    \item{\code{Cynosurus_cristatus}}{}
    \item{\code{Dactylis_glomerata}}{}
    \item{\code{Daucus_carota}}{}
    \item{\code{Festuca_pratensis_sl}}{}
    \item{\code{Geranium_sylvaticum}}{}
    \item{\code{Leontodon_hispidus_sl}}{}
    \item{\code{Potentilla_erecta}}{}
    \item{\code{Pritzelago_alpina_sstr}}{}
    \item{\code{Prunella_vulgaris}}{}
    \item{\code{Ranunculus_acris_sl}}{}
    \item{\code{Saxifraga_oppositifolia}}{}
    \item{\code{Soldanella_alpina}}{}
    \item{\code{Taraxacum_officinale_aggr}}{}
    \item{\code{Trifolium_repens_sstr}}{}
    \item{\code{Veronica_chamaedrys}}{}
    \item{\code{Parnassia_palustris}}{}
    \item{\code{glm_Agrostis_capillaris}}{GLM model for the species Agrostis_capillaris.}
    \item{\code{glm_Leontodon_hispidus_sl}}{GLM model for the species Leontodon_hispidus_sl.}
    \item{\code{glm_Dactylis_glomerata}}{GLM model for the species Dactylis_glomerata.}
    \item{\code{glm_Trifolium_repens_sstr}}{GLM model for the species Trifolium_repens_sstr.}
    \item{\code{glm_Geranium_sylvaticum}}{GLM model for the species Geranium_sylvaticum.}
    \item{\code{glm_Ranunculus_acris_sl}}{GLM model for the species Ranunculus_acris_sl.}
    \item{\code{glm_Prunella_vulgaris}}{GLM model for the species Prunella_vulgaris.}
    \item{\code{glm_Veronica_chamaedrys}}{GLM model for the species Veronica_chamaedrys.}
    \item{\code{glm_Taraxacum_officinale_aggr}}{GLM model for the species Taraxacum_officinale_aggr.}
    \item{\code{glm_Plantago_lanceolata}}{GLM model for the species Plantago_lanceolata.}
    \item{\code{glm_Potentilla_erecta}}{GLM model for the species Potentilla_erecta.}
    \item{\code{glm_Carex_sempervirens}}{GLM model for the species Carex_sempervirens.}
    \item{\code{glm_Soldanella_alpina}}{GLM model for the species Soldanella_alpina.}
    \item{\code{glm_Cynosurus_cristatus}}{GLM model for the species Cynosurus_cristatus.}
    \item{\code{glm_Campanula_scheuchzeri}}{GLM model for the species Campanula_scheuchzeri.}
    \item{\code{glm_Festuca_pratensis_sl}}{GLM model for the species Festuca_pratensis_sl.}
    \item{\code{gbm_Bromus_erectus_sstr}}{GBM model for the species Bromus_erectus_sstr.}
    \item{\code{glm_Saxifraga_oppositifolia}}{GLM model for the species Saxifraga_oppositifolia.}
    \item{\code{glm_Daucus_carota}}{GLM model for the species Daucus_carota.}
    \item{\code{glm_Pritzelago_alpina_sstr}}{GLM model for the species Pritzelago_alpina_sstr.}
    \item{\code{glm_Bromus_erectus_sstr}}{GLM model for the species Bromus_erectus_sstr.}
    \item{\code{gbm_Saxifraga_oppositifolia}}{GBM model for the species Saxifraga_oppositifolia.}
    \item{\code{gbm_Daucus_carota}}{GBM model for the species Daucus_carota.}
    \item{\code{gbm_Pritzelago_alpina_sstr}}{GBM model for the species Pritzelago_alpina_sstr.}
}
}
\details{
The study area is the Western Swiss Alps of Canton de Vaud, Switzerland.

Five topo-climatic explanatory variables to calibrate the SDMs: growing degree days (with a 0 degrees Celsius threshold); moisture index over the growing season (average values for June to August in mm day-1); slope (in degrees); topographic position (an integrated and unitless measure of topographic exposure; Zimmermann et al., 2007); and the annual sum of radiation (in kJ m-2 year-1). The spatial resolution of the predictor is 25 m x 25 m so that the models could capture most of the small-scale variations of the climatic factors in the mountainous areas.

Two modelling techniques were used to produce the SDMs: generalized linear models (GLM; McCullagh & Nelder, 1989; R library 'glm') and generalized boosted models (GBM; Friedman, 2001; R library 'gbm'). The SDMs correpond to 20 species: Agrostis_capillaris, Leontodon_hispidus_sl, Dactylis_glomerata, Trifolium_repens_sstr, Geranium_sylvaticum, Ranunculus_acris_sl, Prunella_vulgaris, Veronica_chamaedrys, Taraxacum_officinale_aggr, Plantago_lanceolata, Potentilla_erecta, Carex_sempervirens, Soldanella_alpina, Cynosurus_cristatus, Campanula_scheuchzeri, Festuca_pratensis_sl, Daucus_carota, Pritzelago_alpina_sstr, Bromus_erectus_sstr and Saxifraga_oppositifolia.}


\references{Guisan, A. 1997. Distribution de taxons vegetaux dans un environnement alpin: Application de modelisations statistiques dans un systeme d'information geographique. PhD Thesis, University of Geneva, Switzerland.

Guisan, A., J.P. Theurillat. and F. Kienast. 1998. Predicting the potential distribution of plant species in an alpine environment. \emph{Journal of Vegetation Science}, \bold{9}, 65-74.

Guisan, A. and J.P. Theurillat. 2000. Assessing alpine plant vulnerability to climate change: A modeling perspective. \emph{Integrated Assessment}, \bold{1}, 307-320.

Guisan, A. and J.P. Theurillat. 2000. Equilibrium modeling of alpine plant distribution and climate change : How far can we go? \emph{Phytocoenologia}, \bold{30}(3-4), 353-384.

Dubuis A., J. Pottier, V. Rion, L. Pellissier, J.P. Theurillat and A. Guisan. 2011. Predicting spatial patterns of plant species richness: A comparison of direct macroecological and species stacking approaches. \emph{Diversity and Distributions}, \bold{17}, 1122-1131.

Zimmermann, N.E., T.C. Edwards, G.G Moisen, T.S. Frescino and J.A. Blackard. 2007. Remote sensing-based predictors improve distribution models of rare, early successional and broadleaf tree species in Utah. \emph{Journal of Applied Ecology} \bold{44}, 1057-1067.
}

\author{Antoine Guisan \email{antoine.guisan@unil.ch}, Anne Dubuis \email{anne.dubuis@gmail.com} and Valeria Di Cola \email{valeria.dicola@unil.ch}
}

\examples{
data(ecospat.testData)
str(ecospat.testData)
dim(ecospat.testData)
names(ecospat.testData)
}
