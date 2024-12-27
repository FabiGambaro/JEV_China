setwd("~/Dropbox/JEV_study")

require(seqinr)
require(plyr)
require(readr)
require(treeio)
require(ape)
require(seraphim)
require(diagram)
require(fields)
require(lubridate)
require(maptools)
require(ncdf4)
require(raster)
require(sp)
require(ggtree)
require(randomcoloR)


# 1. Identifying the different lineages (clades) following introduction event into study area

tree = readAnnotatedNexus("BEAST_dta/JEV_DTA.MCC.tree")
data = read.csv("data/JEV_all_sequences.csv", sep=";", head=T)

CHN_introductions = c() # tree$edge ID corresponding to an introduction event
for (i in 1:dim(tree$edge)[1]){
  if (is.null(tree$annotations[[i]])){
    print(c(h,i))
  }	else		{
    if (tree$annotations[[i]]$location == "China"){
      index = which(tree$edge[,2]==tree$edge[i,1])
      if (tree$annotations[[index]]$location != "China"){
        CHN_introductions = c(CHN_introductions, i)
      }
    }
  }
}
for (i in 1:length(CHN_introductions)){
  if (i == 1) clusters1 = list() # list of sequence IDs related to each introduction event
  if (tree$edge[CHN_introductions[i],2]%in%tree$edge[,1]){
    subtree = treeio::tree_subset(tree, tree$edge[CHN_introductions[i],2], levels_back=0)
    clusters1[[i]] = gsub("'","",subtree$tip.label)
  }	else		{
    clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[CHN_introductions[i],2]])
  }
}
# step to remove the nested clades correspondingb to distinct events
for (i in 2:length(clusters1)) {
  for (j in 1:(i-1)){
    if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]])){
      clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
    }
    if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]])){
      clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
    }
  }
}

for (i in 1:length(CHN_introductions)){
  tab = c()
  if (i == 1) clusters2 = list()
  for (j in 1:length(clusters1[[i]])){
    index = which(data[,"name"]==clusters1[[i]][j])
    if ((length(index) == 1)&&(data[index,"country"]=="mainland China")){
      collection_date = gsub("\\/","-",data[index,"date"])
      prefecture_or_city = data[index,"prefecture.city"]
      if (is.null(prefecture_or_city)){
        prefecture_or_city <- NA
      }
      county = data[index,"county"]
      if (is.null(county)){
        county <- NA
      }
      line = cbind(collection_date, prefecture_or_city, county)
      row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
    }	else	{
       print(clusters1[[i]][j])
       print(data[index,"country"])
    }
  }
  if (!is.null(tab)){
    colnames(tab) = c("collection_date", "prefecture_or_city", "county")
    clusters2[[i]] = tab
  }
}
saveRDS(clusters1, "BEAST_dta/clusters1.rds"); saveRDS(clusters2, "BEAST_dta/clusters2.rds")
clusters1 = readRDS("BEAST_dta/clusters1.rds"); clusters2 = readRDS("BEAST_dta/clusters2.rds")
cluster_sizes = c()
for (i in 1:length(clusters2)){
  cluster_sizes = c(cluster_sizes, dim(clusters2[[i]])[1])
}
cat("\t\tDetected JEV clades circualting in Mainland China :",length(cluster_sizes),"\n",sep="")


## Creating final clusters file
clusters3 <- list()
for (i in 1:length(clusters2)){
  name_clade <- paste0("Clade_", i)
  df <- as.data.frame(clusters2[[i]])
  df$strain <- rownames(df)
  clusters3[[name_clade]] <- df
}
saveRDS(clusters3,"BEAST_dta/clusters3.rds")

final_clusters <- list()
for (i in 1:length(clusters3)){
  name_clade <- names(clusters3)[[i]]
  df <- as.data.frame(clusters3[[i]])
  if (nrow(df) >= 3) {
    final_clusters[[name_clade]] <- df
  } 
}

all_clades <- ldply(final_clusters, rbind)
write_tsv(all_clades, "BEAST_dta/clade_sequences.tsv")
clades_table <- as.data.frame(table(all_clades$.id))
write_tsv(clades_table, "BEAST_dta/table_clades.tsv")
saveRDS(final_clusters,"BEAST_dta/final_clusters.rds")


#2. Plot and color MCC tree showing the different clades the different clades (fig S1)

source("Tree_data_extraction1.r") # for the MCC tree with tip labels
source("Tree_data_extraction2.r") # for posterior trees
source("mccExtractions.r") # for the MCC no tiplabels
source("DTA_tree_extraction1.r")

tree <- read.beast("BEAST_dta/JEV_DTA.MCC.tree")
myplot <- ggtree(tree, mrsd="2020-09-15") + theme_tree2()

clades_sequences <- read_tsv("BEAST_dta/clade_sequences.tsv")
df <- as.data.frame(matrix(nrow = length(tree@phylo$tip.label), ncol = 2))
colnames(df) <- c("tip.label", "clade") 
df$tip.label <- tree@phylo$tip.label
df$clade <- 'no.clade'

matches <- match(df$tip.label, clades_sequences$strain)
df$clade[!is.na(matches)] <- clades_sequences$.id[matches[!is.na(matches)]]
df$clade <- as.factor(df$clade)

df$clade <- factor(df$clade, levels = c("Clade_6", "Clade_9","Clade_13", "Clade_17",
                                        "Clade_19", "Clade_21", "Clade_22", "no.clade"))
mycolors <- c("#E76F51","#F4A261", "#EFB366","#E9C46A","#8AB17D","#2A9D8F","#287271", "#494949")

pdf("BEAST_dta/MCC-clades.pdf")
myplot2 <- myplot %<+% df +
  geom_tippoint(pch=16, aes(color=clade), size=2.5,alpha=.6) + scale_color_manual(values=mycolors)
myplot2
dev.off()



# 3. Running continuous phylogeographic model using 1000 emperical trees from skygrid analysis

## get 1000 emperical trees
log <- read.table("BEAST_Clade_19/Clade_19_skygrid.log", head = T)
burnIn <- round(0.1 * dim(log)[1]) +1
nberOfTreesToSample <- 1000
trees <- scan("BEAST_Clade_19/Clade_19_skygrid.trees",what="", sep="\n", quiet=T, blank.lines.skip=F)

index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]
index2 = index1 + burnIn + 1
indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
selected_trees = c(trees[c(1:index1,indices)],"End;")
write(selected_trees,"Clade_19_skygrid_1000.trees")


#Add polygons as homogeneous prior of sampling locations

data <- read_tsv("BEAST_Clade_19/Clade19_data.tsv")
directoryPolygons = "Sample_polygons"
xml = scan("BEAST_Clade_19/Clade_19_cauchy_template_polygons_clean.xml", what="list", sep="\n", quiet=T, blank.lines.skip=F)
sink(file = "BEAST_Clade_19/Clade_19_cauchy_polygons_wgs84_clean.xml")

for (h in 1:length(xml)){
  cat(xml[h]); cat("\n")
  
  if (xml[h] == "\t<!-- INSERT leafTraitParameter -->"){
    for (i in 1:nrow(data)){
      sequence_name <- data$name[i]
      cat(paste0("\t<leafTraitParameter id=\"", sequence_name,".trait\" taxon=\"",
                 sequence_name,"\">")) ;cat("\n")
      cat(paste0("\t\t<treeModel idref=\"", "treeModel","\"/>")) ;cat("\n")
      cat(paste0("\t\t<parameter idref=\"", "leaf.location", "\"/>"))
      cat("\n")
      cat(paste0("\t</leafTraitParameter>")); cat("\n")
    }
    cat("\n")
    for (i in 1:nrow(data)){
      sequence_name <- data$name[i]
      cat(paste0("\t<flatGeoSpatialPrior id=\"",sequence_name, "_polygon\" taxon=\"", 
                 sequence_name,  "\" kmlFileName=\"", directoryPolygons,"/", sequence_name,
                 ".kml\" inside=\"true\" union=\"true\" cache=\"true\">"))
      cat("\n")
      cat(paste0("\t\t<data>")); cat("\n")
      cat(paste0("\t\t\t<parameter idref=\"",sequence_name,".trait\"/>")); cat("\n")
      cat(paste0("\t\t</data>")); cat("\n")
      cat(paste0("\t</flatGeoSpatialPrior>")); cat("\n")
    }
    cat("\n")
  } 
  if (xml[h]=="\t\t<!-- INSERT uniformGeoSpatialOperator -->")
  {
    cat("\n")
    for (i in 1:nrow(data)){
      sequence_name <- data$name[i]
      cat(paste0("\t\t<uniformGeoSpatialOperator weight=\"0.01\">")); cat("\n")
      cat(paste0("\t\t\t<parameter idref=\"",sequence_name,".trait\"/>"));cat("\n")
      cat(paste0("\t\t\t<flatGeoSpatialPrior idref=\"",sequence_name,"_polygon\"/>"))
      cat("\n")
      cat(paste0("\t\t</uniformGeoSpatialOperator>")); cat("\n")
    }
    cat("\n")
  }
  if (xml[h]=="\t\t\t\t<!-- INSERT geoDistributionCollection -->")
  {
    cat("\n")
    cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
    for (i in 1:nrow(data)){
      sequence_name <- data$name[i]
      cat(paste0("\t\t\t<flatGeoSpatialPrior idref=\"",sequence_name,"_polygon\"/>"))
      cat("\n")
    }
    cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
    cat("\n")
  }
  
}

sink(NULL)


#4. Analysis of the continous phylogeographic analysis

## Building the MCC tree 
log <- read_delim("BEAST_Clade_19/Clade_19_cauchy_polygons_wgs84.log", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 3)
burnIn <- round(0.1 * dim(log)[1]) +1
system(paste0("/Applications/BEAST_1/bin/treeannotator -burninTrees ",burnIn," -heights keep ",
              "BEAST_Clade_19/Clade_19_cauchy_polygons_wgs84.trees", " ", "BEAST_Clade_19/Clade_19_cauchy_polygons_wgs84.MCC.tree"))



## Extracting spatio-temporal information embedded in MCC and posterior trees of the RRW analysis

Directory = "BEAST_Clade_19/"
nberOfExtractionFiles <- 900
nberOfTreesToSample = nberOfExtractionFiles; randomSampling = FALSE; coordinateAttributeName = "location"; nberOfCores = 10
treeFiles <- list.files("BEAST_Clade_19/",pattern = ".trees")
treeFiles <- treeFiles[1]

for (j in 1:length(treeFiles)){
  tree_name <- gsub(".trees","",treeFiles[j])
  mostRecentSamplingDatum <- 2020.6666666666667 # 2020-09
  mcc_tre <- readAnnotatedNexus(paste0(Directory,tree_name,"_MCC.tree"))
  mcc_tab = as.data.frame(mccExtractions(mcc_tre, mostRecentSamplingDatum))
  write.csv(mcc_tab, paste0(Directory, tree_name,".csv"), row.names=F, quote=F)
  
  allTrees <- scan(file=paste0(Directory,tree_name,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
  localTreesDirectory <- paste0(Directory,tree_name, "_ext")
  treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
  
}


##estimating the HPD region for each time slice
startDatum = min(mcc_tab[,"startYear"])
percentage = 80; prob = percentage/100
precision = 5 #units years

nberOfExtractionFiles <- 900
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

admin2 = shapefile("Chinese_shapefiles/Admin-2_polygons.shp")

#defining the different colour scales to use
colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[26:126])
minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])

startYears_indices = (((mcc_tab[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
startYears_colours = colourScale[startYears_indices]
endYears_colours = colourScale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons)){
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colourScale[polygon_index],"40") 
}

pdf(paste0(Directory,"RRW_dispersal_plot/",tree_name,".pdf"))

admins2 = shapefile("Chinese_shapefiles/Admin-2_polygons.shp")
admins2 = spTransform(admins2, CRS("+init=epsg:4326")); admins1 = rgeos::gUnaryUnion(admins2)
e_nodes = extent(raster("Environmental_files/Elevation_JEV_04.tif"))
admins1_cropped = crop(admins1, e_nodes); admins2_cropped = crop(admins2, e_nodes)

plot(admins2_cropped, col="#6C757D", border=NA, lwd=0.01)
plot(admins2_cropped, add=T, lwd=0.01, border="white")
#text(admin3, labels = admin3@data$NAME, col = "black", cex = 0.6)

for (i in 1:length(polygons)) {
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}

for (i in 1:dim(mcc_tab)[1]) {
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]),
              cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA,
              arr.pos=F, curve=0.1, dr=NA, endhead=F)
}

for (i in 1:dim(mcc_tab)[1]){
  if (i == 1){
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16, col=colourScale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1, col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8)
}

rect(xmin(admins2_cropped), ymin(admins2_cropped), xmax(admins2_cropped),
     ymax(admins2_cropped), xpd=T, lwd=0.2)
axis(1, c(ceiling(xmin(admins2_cropped)), floor(xmax(admins2_cropped))),
     pos=ymin(admins2_cropped), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(floor(ymin(admins2_cropped)), floor(ymax(admins2_cropped))),
     pos=xmin(admins2_cropped), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=1, tck=-0.01, col.axis="gray30")

selectedDates = decimal_date(ymd(c("1980-01-01","1990-01-01","2000-01-01","2010-01-01", "2020-01-01")))
selectedLabels =c("1980","1990","2000","2010","2020")

rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
plot(rast, legend.only = T, add = T, col = colourScale, legend.width = 0.5, legend.shrink = 0.3,
     smallplot = c(0.40, 0.80, 0.14, 0.155), #dimensions of the legend plot
     legend.args = list(text = "", cex = 0.7, line = 0.3, col = "gray30"),
     horizontal = T, #legend horizontal
     axis.args = list(cex.axis = 0.6, lwd = 0, lwd.tick = 0.2, tck = -0.5, col.axis = "gray30",line = 0, mgp = c(0, -0.02, 0), at=selectedDates, labels=selectedLabels))
dev.off()



pdf("temporal_spread.pdf", width=10, height=5.8) # fig S2
par(mfrow=c(2,4), mar=c(0,0,0,0), oma=c(1,1,1,1), mgp=c(0,0.4,0), lwd=0.2, bty="o")

cutOffs = c(1980, 2000, 2015, maxYear); croppingPolygons = FALSE
for (h in 1:length(cutOffs)){
  plot(admins1_cropped, col="gray95", border=NA, lwd=0.01, axes=F)
  rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc_tab[,"endYear"])
  if (h == length(cutOffs)){
    plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.150,0.510,0.900,0.915),
         legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
         axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.9, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.10,0)))
  }
  for (i in 1:length(polygons)){
    if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
    {
      for (j in 1:length(polygons[[i]]@polygons))
      {
        polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
      }
      pol = polygons[[i]]; crs(pol) = crs(admins1)
      if (croppingPolygons == TRUE) pol = crop(pol, admins1)
      plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
    }
  }
  plot(admins2_cropped, col=NA, border="gray75", lwd=0.2, add=T)
  plot(admins1_cropped, col=NA, border="gray50", lwd=0.4, add=T)
  for (i in 1:dim(mcc_tab)[1])
  {
    if (mcc_tab[i,"endYear"] <= cutOffs[h])
    {
      curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
                  arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
    }
  }
  for (i in dim(mcc_tab)[1]:1)
  {
    if ((mcc_tab[i,"startYear"] <= cutOffs[h])&(!mcc_tab[i,"node1"]%in%mcc_tab[,"node2"]))
    {
      startYears_index = (((mcc_tab[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
      points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16, col=colourScale[startYears_index], cex=0.7)
      points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
    }
    if (mcc_tab[i,"endYear"] <= cutOffs[h])
    {
      points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
      points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
    }
  }
  if (h < length(cutOffs)) mtext(cutOffs[h], side=1, adj=0, line=-18, at=-300000, cex=0.7, font=1, col="gray30")
}
dev.off()


# 9. Dispersal statistics based on the continuous phylogeographic reconstructions

tree_name <- gsub(".trees","",treeFiles[j])
localTreesDirectory <- paste0(Directory,tree_name, "_ext")

timeSlices = 100; onlyTipBranches = FALSE; showingPlots = TRUE; nberOfCores = 5; slidingWindow = 1
outputName = paste0(Directory,"RRW_dispersal_statistics/", tree_name)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)

stats = read.table(paste0(Directory,"RRW_dispersal_statistics/", tree_name, "_estimated_dispersal_statistics.txt"), head=T) #load the files 
wldv = stats[,"weighted_branch_dispersal_velocity"]; wldv_md = round(median(wldv),2)
wldv_qs = round(quantile(wldv,c(0.025,0.975)),2); wldv_hpd = round(HDInterval::hdi(wldv)[1:2],2)
wdc = stats[,"weighted_diffusion_coefficient"]; wdc_md = round(median(wdc),1)
wdc_qs = round(quantile(wdc,c(0.025,0.975)),1); wdc_hpd = round(HDInterval::hdi(wdc)[1:2],2)
cat(tree_name,":\tWLDV = ",wldv_md," km/day [",wldv_hpd[1],"-",wldv_hpd[2],"]\t\tWDC = ",wdc_md," km2/year [",wdc_hpd[1],"-",wdc_hpd[2],"]","\n",sep="")	

df_wldv <- as.data.frame(stats[,"weighted_branch_dispersal_velocity"])
colnames(df_wldv) <- "wldv"

ggplot(df_wldv, aes(x = wldv)) +
  geom_density(color = "black")

df_wdc <- as.data.frame(stats[,"weighted_diffusion_coefficient"])
colnames(df_wdc) <- "wdc"

pdf(paste0(Directory,"RRW_dispersal_statistics/", tree_name,"_wdc_plot.pdf"))
ggplot(df_wdc, aes(x = wdc)) +
  geom_density(color = "black", fill = "orange", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 40000))+
  xlab("wdc km2/day")+
  theme_bw() 
dev.off()


# 10. PLotting the skygrid and coloring according the colors used in the continous RRW model. 

log <- read_delim("BEAST_Clade_19/Clade_19_skygrid_per_month.log", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 4)
burnin_fraction <- 0.1
burnIn <- round(burnin_fraction * dim(log)[1]) +1
log <- log[,c(1,9:257)] #keep only pop.size values
log1 = log[burnIn:dim(log)[1],]
selected_cols <- grep("^skygrid.logPopSize", colnames(log1), value = TRUE)

names_vector <- c()
mean_vector <- c()
low_hpd_vector <- c()
high_hpd_vector <- c()

for (col_name in selected_cols) {
  mean_val <- mean(log1[[col_name]])
  hpd_low <- HDInterval::hdi(log1[[col_name]])[1]
  hpd_high <-  HDInterval::hdi(log1[[col_name]])[2]
  #store results
  names_vector <- c(names_vector, col_name)
  mean_vector <- c(mean_vector, mean_val)
  low_hpd_vector <- c(low_hpd_vector, hpd_low)
  high_hpd_vector <- c(high_hpd_vector,hpd_high)
}

result_df <- data.frame(
  name = names_vector,
  mean = mean_vector,
  low_HPD = low_hpd_vector,
  high_HPD = high_hpd_vector
)

#get grid points that I used to build the skygrid (grid.point.R script)
separator <- "XXXX"
grid.points.Jan <- paste0(2000:2020, c("-01-01"))
grid.points.Feb <- paste0(2000:2020, c("-02-01"))
grid.points.Mar <- paste0(2000:2020, c("-03-01"))
grid.points.Apr <- paste0(2000:2020, c("-04-01"))
grid.points.May <- paste0(2000:2020, c("-05-01"))
grid.points.Jun <- paste0(2000:2020, c("-06-01"))
grid.points.Jul <- paste0(2000:2020, c("-07-01"))
grid.points.Aug <- paste0(2000:2020, c("-08-01"))
grid.points.Sep <- paste0(2000:2020, c("-09-30"))
grid.points.Oct <- paste0(2000:2020, c("-10-01"))
grid.points.Nov <- paste0(2000:2020, c("-11-01"))
grid.points.Dec <- paste0(2000:2020, c("-12-01"))

combo <- paste(grid.points.Jan, grid.points.Feb,grid.points.Mar, grid.points.Apr,
               grid.points.May,grid.points.Jun, grid.points.Jul,grid.points.Aug,
               grid.points.Sep,grid.points.Oct,grid.points.Nov,grid.points.Dec,
               sep = separator)

grid.points <- unlist(strsplit(x = combo, split = separator))
grid.points <- grid.points[-c(250:252)]
grid.points <- c(head(grid.points, -1))
grid.points.decimal <- lubridate::decimal_date(as.Date(grid.points))
MRSD <- lubridate::decimal_date(as.Date("2020-09-15")) 
length( grid.points.decimal)
grid.points.decimal2 <- c(grid.points.decimal, MRSD)
grid.points.decimal2 <- rev(grid.points.decimal2)

df <- cbind(result_df, grid.points.decimal2)
colnames(df) <- c("name","mean","low_HPD","high_HPD","time")


pdf(paste0(Directory,"Clade_19_per_month_skygrid.pdf"), width=10.0, height=1.7)
par(mar=c(1.7,2.7,0.7,1), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o")
xMin = min(df[,"time"]); xMax = max(df[,"time"]); xMin = minYear; 
timeSlice = df[1,"time"]-df[2,"time"]
yMin = min(df[,"low_HPD"]); yMax = max(df[,"high_HPD"])
yMin = 0; yMax = 8 # visualisation parameters for the plot (to be edited)
colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[26:126])
colours = colourScale[(((df[,c("time")]-minYear)/(maxYear-minYear))*100)+1]

plot(df[,"time"], df[,"mean"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
xx_l = c(df[,c("time")],rev(df[,c("time")]));  yy_l = c(df[,"low_HPD"],rev(df[,"high_HPD"]))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)

for (i in 1:length(df[,"time"])){
  x1 = df[i,"time"]-(timeSlice/2); x2 = df[i,"time"]+(timeSlice/2)
  y1 = df[i,"low_HPD"]-0.2; y2 = df[i,"high_HPD"]+0.2
  polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colours[i], border=NA)
}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
lines(df[,"time"], df[,"mean"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col="gray30", col.axis="gray30", mgp=c(0,0.04,0), at=seq(2000,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col="gray30", col.axis="gray30", mgp=c(1,0.30,0), at=seq(0,8,1))
mtext("log(Ne)", side=2, col="gray30", cex=0.50, line=1.2, las=3)
dev.off()


# 11. Plotting the different environmental rasters

admins2 = shapefile("Chinese_shapefiles/Admin-2_polygons.shp")
admins2 = spTransform(admins2, CRS("+init=epsg:4326")); admins1 = rgeos::gUnaryUnion(admins2)
e_nodes = extent(raster("Environmental_files/Elevation_JEV_08.tif"))
admins1_cropped = crop(admins1, e_nodes); admins2_cropped = crop(admins2, e_nodes)

envVariableFiles = c("Annual_mean_temperature","Annual_precipitation","ENM_C_tritaeniorhynchus","Pig_pop_density_log",
                     "Land_cover_water", "total_rice", "Human_pop_density_log","Land_cover_urban_areas")

envVariableNames1 = c("Annual","Annual","Culex","Pig pop.","","","Human pop","")
envVariableNames2 = c("mean temp.","precipitation","tritaeniorhynchus","density (log)", 
                      "Water bodies", "Rice fields", "density (log)","Urban areas")


#load
rS = list(); cols = list(); colour1 = "gray98"
for (i in 1:length(envVariableFiles)){
  rS[[i]] <- raster(paste0("Environmental_files/", envVariableFiles[i], "_JEV_08.tif"))
  crs(rS[[i]]) = CRS("+init=epsg:4326"); rS[[i]] = crop(mask(crop(projectRaster(rS[[i]],crs=crs(admins2)),admins1),admins1),e_nodes)
}

#urban areas
rS[[8]][!is.na(rS[[8]][])] = log(rS[[8]][!is.na(rS[[8]][])]+1) # log-transformation of urban areas
#precipitation
rS[[2]][] = rS[[2]][]/100 # legend: temperature in °C and precipitation in meters
#urban areas
cols[[8]] = c("#FAFAFA",colorRampPalette(brewer.pal(9,"BuPu"))(100))
#temp
cols[[1]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
#precipitation
cols[[2]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
#culex
cols[[3]] = colorRampPalette(brewer.pal(9,"PuRd"))(100)
#pig
cols[[4]] = colorRampPalette(brewer.pal(9,"RdPu"))(120)[1:100]
#water
colour2 = "deepskyblue4"; r = rS[[5]]
cols[[5]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
#total rice
cols[[6]] = c("#FAFAFA",colorRampPalette(brewer.pal(9,"Greens"))(99))
#human pop size
cols[[7]] = colorRampPalette(brewer.pal(9,"PuBu"))(100)


pdf("Environmental_files/selected_factors.pdf", width=10, height=4.2); colNA = "white"
par(mfrow=c(2,4), mar=c(0,0,0,0), oma=c(1,1,1.5,1), mgp=c(0,0.4,0), lwd=0.2, bty="o", col="gray30")
for (i in 1:length(rS)){
  plot(admins1_cropped, col=NA, border=NA, lwd=0.001)
  plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA=colNA, add=T)
  plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.830,0.845,0.100,0.500),
       legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=F,
       axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.50,0)))
  if (nchar(envVariableNames1[i] > 0)) mtext(envVariableNames1[i], side=1, adj=0, line=-14, at=93, cex=0.7, font=1, col="gray30")
  if (nchar(envVariableNames2[i] > 0)) mtext(envVariableNames2[i], side=1, adj=0, line=-13, at=93, cex=0.7, font=1, col="gray30")
  plot(admins2_cropped, col=NA, border="gray75", lwd=0.2, add=T)
  plot(admins1_cropped, col=NA, border="gray50", lwd=0.4, add=T)
}
dev.off()


# 12. Generating a null model of lineages dispersal

nberOfExtractionFiles = 900
localTreesDirectory <- paste0(Directory,tree_name, "_ext")
if (!file.exists("Convex_hull_rast.tif"))
{
  points1 = c(); points2 = c()
  for (i in 1:nberOfExtractionFiles)
  {
    tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
    points1 = rbind(points1, tab[,c("startLon","startLat")])
    points2 = rbind(points2, tab[,c("endLon","endLat")])
  }
  colnames(points1) = c("lon","lat"); colnames(points2) = c("lon","lat")
  points = rbind(points1, points2); jitter = 0
  if (jitter > 0)
  {
    points1 = points; points2 = points; points3 = points; points4 = points
    points1[,1] = points1[,1]-jitter; points2[,1] = points2[,1]+jitter
    points3[,2] = points3[,2]-jitter; points4[,2] = points4[,2]+jitter
    points = rbind(points, points1, points2, points3, points4)
  }
  hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
  ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
  rast1 = raster("Environmental_files/Elevation_JEV_04.tif")
  crs(rast1) = CRS("+init=epsg:4326")
  rast2 = projectRaster(rast1, crs=crs(admins2))
  rast3 = rast2; rast3[!is.na(rast3[])] = 0
  rast4 = rasterize(sps, rast3, getCover=T, fun=avg)
  rast5 = rast4; rast5[rast5[]==0] = NA
  rast5[is.na(rast3[])] = NA; backgroundRaster = rast5
  backgroundRaster = crop(mask(rast5,admins2),admins2)
  writeRaster(backgroundRaster, "Convex_hull_rast.tif", overwrite=T)
}	else		{
  backgroundRaster = raster("Convex_hull_rast.tif")
}
envVariables = list(backgroundRaster); randomProcedure = 3; nberOfCores = 5
treesRandomisation(localTreesDirectory, nberOfExtractionFiles, envVariables, randomProcedure, nberOfCores)

library(rgeos)
library(seraphim)

# 13. Testing the impact of environmental factors on lineage dispersal velocity

localTreesDirectory <- paste0(Directory,tree_name, "_ext"); nberOfExtractionFiles = 100
envVariableFiles = c("Land_cover_savannas","Land_cover_forests","Land_cover_croplands","Land_cover_urban_areas",
                     "Annual_mean_temperature","Annual_precipitation","ENM_C_tritaeniorhynchus","Pig_pop_density_log",
                     "Land_cover_wetlands", "Land_cover_water", "total_rice", "irrigated_rice")

nberOfExtractionFiles = 100; nberOfRandomisations = 0; randomProcedure = 3
showingPlots = FALSE; nberOfCores = 10; OS = "Unix"; randomisations = FALSE; fourCells = FALSE; c = 0
envVariables = list(); resistances = list(); avgResistances = list()

#loading my raster files into a list and creating a list which will indicate if the raster will be treated as resistance ir not
for (k in c(10,100,1000)){
  for (i in 1:length(envVariableFiles)){
    c = c+1
    rast = raster(paste("Environmental_files/",envVariableFiles[i],"_JEV_08.tif",sep=""))
    rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
    names(rast) = paste(envVariableFiles[i],"_k",k,sep="")
    envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableFiles[i],"_k",k,sep="")
    resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
  }
  for (i in 1:length(envVariableFiles)){
    c = c+1
    rast = raster(paste("Environmental_files/",envVariableFiles[i],"_JEV_08.tif",sep=""))
    rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
    names(rast) = paste(envVariableFiles[i],"_k",k,sep="")
    envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableFiles[i],"_k",k,sep="")
    resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
  }
}

pathModel = 3; outputName = paste0("JEV_CS_extractions")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
              nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)	

pathModel = 3; outputName = paste0("JEV_CS_randomisations")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
              nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=T)

pathModel = 2; outputName = paste0("JEV_LC_extractions")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
              nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)
pathModel = 2; outputName = paste0("JEV_LC_randomisations")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
              nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=T)

Qs_tab = c(); res = c("R","C")
tab = read.table(paste0("SpreadFactorsOutput/JEV_CS_extractions_linear_regression_results.txt"), head=T)
for (i in 1:length(envVariableFiles)){
  for (j in 1:length(res)){
    for (k in c(10,100,1000)){
      Qs = tab[,which(colnames(tab)==paste0("Univariate_LR_delta_R2_",envVariableFiles[i],"_k",k,"_",res[j]))]
      if (sum(Qs>0) >= 90) print(paste0("Univariate_LR_delta_R2_",envVariableFiles[i],"_k",k,"_",res[j]))
      Qs_tab = rbind(Qs_tab, cbind(round(sum(Qs>0)/length(Qs),2), paste0("Univariate_LR_delta_R2_",envVariableFiles[i],"_k",k,"_",res[j])))
    }
  }
}
selected_variables = Qs_tab[which(as.numeric(Qs_tab[,1])>0.9),2] #none with CS or LC


# 14. Testing the impact of environmental factors on lineage dispersal locations

localTreesDirectory <- paste0(Directory,tree_name, "_ext"); nberOfExtractionFiles = 100; envVariables = list()
envVariableFiles = c("Land_cover_savannas","Land_cover_forests","Land_cover_croplands","Land_cover_urban_areas",
                     "Annual_mean_temperature","Annual_precipitation","ENM_C_tritaeniorhynchus","Pig_pop_density_log","Human_pop_density_log")

for (i in 1:length(envVariableFiles)){
  envVariables[[i]] = raster(paste("Environmental_files/",envVariableFiles[i],"_JEV_08.tif",sep=""))
}

for (i in 1:nberOfExtractionFiles){
  obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
  ran = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",i,".csv"), header=T)
  envValues_obs = matrix(nrow=dim(obs)[1], ncol=length(envVariables))
  envValues_ran = matrix(nrow=dim(ran)[1], ncol=length(envVariables))
  colnames(envValues_obs) = envVariableFiles
  colnames(envValues_ran) = envVariableFiles
  for (j in 1:length(envVariables)){
    envValues_obs[,j] = raster::extract(envVariables[[j]], SpatialPoints(obs[,c("endLon","endLat")]))
    envValues_ran[,j] = raster::extract(envVariables[[j]], SpatialPoints(ran[,c("endLon","endLat")]))
  }
  write.csv(envValues_obs, paste0(localTreesDirectory,"/EnvValues_obs",i,".csv"), row.names=F, quote=F)
  write.csv(envValues_ran, paste0(localTreesDirectory,"/EnvValues_ran",i,".csv"), row.names=F, quote=F)
}
BFs = matrix(nrow=length(envVariables), ncol=2)
row.names(BFs) = envVariableFiles; colnames(BFs) = c("lower","higher")
meanEnvValues_obs_list = list(); meanEnvValues_ran_list = list()
for (i in 1:length(envVariables)){
  lowerEnvValues = 0; meanEnvValues_obs_list = list(); meanEnvValues_ran_list = list()
  meanEnvValues_obs = rep(NA, nberOfExtractionFiles); meanEnvValues_ran = rep(NA, nberOfExtractionFiles)
  for (j in 1:nberOfExtractionFiles){
    meanEnvValues_obs[j] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_obs",j,".csv"))[,gsub(".tif","",gsub("-",".",envVariableFiles[i]))], na.rm=T)
    meanEnvValues_ran[j] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_ran",j,".csv"))[,gsub(".tif","",gsub("-",".",envVariableFiles[i]))], na.rm=T)
    if (meanEnvValues_obs[j] < meanEnvValues_ran[j]) lowerEnvValues = lowerEnvValues+1				
  }
  p = lowerEnvValues/nberOfExtractionFiles; BFs[i,"lower"] = round((p/(1-p))/(0.5/(1-0.5)),1)
  p = (1-(lowerEnvValues/nberOfExtractionFiles)); BFs[i,"higher"] = round((p/(1-p))/(0.5/(1-0.5)),1)
  meanEnvValues_obs_list[[i]] = meanEnvValues_obs; meanEnvValues_ran_list[[i]] = meanEnvValues_ran
}
write.csv(BFs, paste0("SpreadFactorsOutput/Seraphim_JEV_rice_water.csv"), quote=F)


# Plot short vs long spread events

nberOfExtractionFiles = 900
Trees_distance <- list()
for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  tab$time <- tab$endYear -tab$startYear
  r1 = tab[,c("startLat", "startLon")]
  r2 = tab[,c("endLat", "endLon")]
  tab$dists = fields::rdist.earth.vec(r1, r2, miles = F)
  tab$vel = (tab$dists/tab$time)/365
  Trees_distance[[i]] <- tab
}

combined_df <- do.call(rbind, Trees_distance)
ggplot(combined_df, aes(x = dists, y = vel)) +
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis") +
  xlab("distance (km)") + 
  ylab("velocity (km/day)") +
  theme_bw()

# Plot kernel density estimate 
quantiles <- quantile(combined_df$dists, probs = c(0.5, 0.95)) # Calculate quantiles
density_values <- density(combined_df$dists) # Calculate density values

pdf("BEAST_Clade_19/RRW_dispersal_plot/distance.pdf", width = 15, height = 5)
ggplot(combined_df, aes(x = dists)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values$x, y = density_values$y),
              aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.3, color = NA) +
  labs(x = "Distance in (km)", y = "Density") +
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  #scale_x_continuous(limits = c(0, 400))+
  theme_bw()
dev.off()

# Create separate datasets for each percentile range
density_values_below_50 <- data.frame(x = density_values$x[density_values$x <= quantiles[1]],
                                      y = density_values$y[density_values$x <= quantiles[1]])
density_values_between_50_95 <- data.frame(x = density_values$x[density_values$x > quantiles[1] & density_values$x <= quantiles[2]],
                                           y = density_values$y[density_values$x > quantiles[1] & density_values$x <= quantiles[2]])
density_values_above_95 <- data.frame(x = density_values$x[density_values$x > quantiles[2]],
                                      y = density_values$y[density_values$x > quantiles[2]])

# Plot each density range separately with different colors
pdf("BEAST_Clade_19/RRW_dispersal_plot/distance2.pdf", width = 15, height = 5)
ggplot(combined_df, aes(x = dists)) +
  geom_density(color = "black") +
  geom_ribbon(data = density_values_below_50, aes(x = x, ymin = 0, ymax = y), fill = "red4", alpha = 0.8, color = NA) +
  geom_ribbon(data = density_values_between_50_95, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.8, color = NA) +
  geom_ribbon(data = density_values_above_95, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.3, color = NA) +
  labs(x = "Distance in (km)", y = "Density") +
  theme_bw()

dev.off()
