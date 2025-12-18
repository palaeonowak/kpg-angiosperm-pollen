#######
# Requires the following files in the working directory:
# "new morphology.csv"
# "new ranges + groups.csv"
# "Lupia1999 morphology.csv"
# "Lupia1999 ranges + groups.csv"
# "new occurrences KPg.csv"
#######

# Load packages
library(Claddis)
library(dispRity)
library(ggplot2)
library(gridExtra)
library(maps)


# Define functions

disp_sv <- function(scores_table, groups_list, n_bootstraps){
  set.seed(1234)
  disp <- dispRity(boot.matrix(custom.subsets(data = scores_table,
                                              group = groups_list),
                               bootstraps = n_bootstraps),
                   metric = c(sum, variances))
  return(disp)
}

summarise_disp<-function(disp, batch_title){
  exclude<-sapply(1:length(disp$disparity), \(x) any(is.na(disp$disparity[[x]][[2]])))
  disp$subsets<-disp$subsets[!exclude]
  disp$disparity<-disp$disparity[!exclude]
  df<-summary(disp)
  df$batch<-batch_title
  df<-df[,c(ncol(df), 1:(ncol(df)-1))]
  return(df)
}

calc_survival<-function(group_list, batch_title){
  n_groups<-length(group_list)/2
  pre<-1:n_groups
  surv_df<-data.frame(batch=batch_title,
                      grouping=names(group_list)[pre],
                      Maastrichtian=sapply(pre, \(x) length(group_list[[x]])),
                      Palaeocene=sapply(n_groups+pre, \(x) length(group_list[[x]])),
                      survivors=sapply(pre, \(x) length(intersect(group_list[[x]], group_list[[n_groups+x]]))),
                      extinct=sapply(pre, \(x) sum(!(group_list[[x]] %in% group_list[[n_groups+x]]))),
                      originating=sapply(pre, \(x) sum(!(group_list[[n_groups+x]] %in% group_list[[x]]))),
                      group_fraction=sapply(pre, \(x) length(group_list[[x]]))/nrow(ranges)*100)
  surv_df$survival_rate<-surv_df$survivors/surv_df$Maastrichtian*100
  surv_df$extinction_rate<-surv_df$extinct/surv_df$Maastrichtian*100
  return(surv_df)
}

group_mannwu<-function(disp, batch_title){
  mwu<-data.frame(batch=batch_title,
                  group=names(groups),
                  W=0,
                  p=0)
  for(i in 1:length(groups)){
    if(!any(is.na(disp[[i]][[2]])) & !any(is.na(disp[[i+length(groups)]][[2]]))){
      mannwhitney<-wilcox.test(disp[[i]][[2]], disp[[i+length(groups)]][[2]])
      mwu$W[i]<-mannwhitney$statistic
      mwu$p[i]<-mannwhitney$p.value
    }else{
      mwu$W[i]<-NA
      mwu$p[i]<-NA
    }
  }
  return(mwu)
}

ggplot_div<-function(data, batch, title){
  data<-data[data$batch==batch,]
  ggplot(data = data, aes(x = Group, y = n, fill = Age))+
    geom_col(position = position_dodge(), colour="black")+
    scale_fill_manual(values = c("#88C96F","#F9AA6F"), breaks = c("Maastrichtian","Palaeocene"), aesthetics = "fill") +
    xlab(label="Pollen groups")+
    ylab(label="Taxonomic diversity")+
    theme_bw()+
    theme(axis.ticks=element_blank())+
    labs(title = title)
}

ggplot_survival<-function(data, batch, title){
  data<-data[data$batch==batch,]
  ggplot(data = data, aes(x=Group, y=value, fill = Turnover))+
    geom_col(position = "stack", colour="black")+
    scale_fill_manual(values = c("#C42503FF" ,"#32AEF6FF"), breaks = c("Extinction", "Survival"), aesthetics = "fill") +
    geom_hline(yintercept = data$value[data$Group=="All" & data$Turnover=="Survival"], linetype=2)+
    xlab(label="Pollen groups")+
    ylab(label="%")+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    theme_bw()+
    theme(axis.ticks=element_blank())+
    labs(title = title)
}

plot_survivors<-function(pcoa_scores, pcoa_varexp, axis_1, axis_2, survivors_id, title){
  plot(pcoa_scores[, axis_1], pcoa_scores[, axis_2],
       type = "n", las = 1, mgp = c(2.5, 0.7, 0),
       xlab = paste("PCoA ", axis_1, " (", round(pcoa_varexp[axis_1], 2), "% variance explained)", sep = ""),
       ylab = paste("PCoA ", axis_2, " (", round(pcoa_varexp[axis_2], 2), "% variance explained)", sep = ""),
       main = title, xaxt = "n", yaxt = "n")
  
  points(pcoa_scores[survivors_id$Maastrichtian, axis_1],
         pcoa_scores[survivors_id$Maastrichtian, axis_2],
         pch = 19, col = "grey20", cex = 1.5)
  points(pcoa_scores[survivors_id$Boundary, axis_1],
         pcoa_scores[survivors_id$Boundary, axis_2],
         pch = 19, col = "green", cex = 1.5)
  points(pcoa_scores[survivors_id$Palaeocene, axis_1],
         pcoa_scores[survivors_id$Palaeocene, axis_2],
         pch = 19, col = "#F9AA6F", cex = 1.5)
  
  axis(side = 1, lwd = 0, lwd.ticks = 1, tcl = 0.3, mgp = c(2.5, 0.5, 0))
  axis(side = 2, lwd = 0, lwd.ticks = 1, tcl = 0.3, mgp = c(2.5, 0.5, 0), las = 1)
}

plot_pcoa_groups<-function(pcoa_scores, pcoa_varexp, axis_1, axis_2, sel_groups, title){
  plot(pcoa_scores[, axis_1], pcoa_scores[, axis_2],
       type = "n", las = 1, mgp = c(2.5, 0.7, 0),
       xlab = paste("PCoA ", axis_1, " (", round(pcoa_varexp[axis_1], 2), "% variance explained)", sep = ""),
       ylab = paste("PCoA ", axis_2, " (", round(pcoa_varexp[axis_2], 2), "% variance explained)", sep = ""),
       main = title, xaxt = "n", yaxt = "n")
  
  points(pcoa_scores[sel_groups$Monosulcate, axis_1],
         pcoa_scores[sel_groups$Monosulcate, axis_2],
         pch = 15, col = "#A2FC3CFF", cex = 1.5)
  
  points(pcoa_scores[sel_groups$Triprojectate, axis_1],
         pcoa_scores[sel_groups$Triprojectate, axis_2],
         pch = 18, col = "#455BCDFF", cex = 1.5)
  
  points(pcoa_scores[setdiff(sel_groups$Colpate, sel_groups$Triprojectate), axis_1],
         pcoa_scores[setdiff(sel_groups$Colpate, sel_groups$Triprojectate), axis_2],
         pch = 18, col = "grey70", cex = 1.5)
  
  points(pcoa_scores[sel_groups$Porate, axis_1],
         pcoa_scores[sel_groups$Porate, axis_2],
         pch = 17, col = "#C42503FF", cex = 1.5)
  
  points(pcoa_scores[sel_groups$Other, axis_1],
         pcoa_scores[sel_groups$Other, axis_2],
         pch = 19, col = "grey90", cex = 1.5)
  
  axis(side = 1, lwd = 0, lwd.ticks = 1, tcl = 0.3, mgp = c(2.5, 0.5, 0))
  axis(side = 2, lwd = 0, lwd.ticks = 1, tcl = 0.3, mgp = c(2.5, 0.5, 0), las = 1)
}

ggplot_disp<-function(disparity, title){
  disparity<-disparity[!sapply(1:length(disparity), \(x) any(is.na(disparity[[x]][[2]])))]
  plotmat<-as.data.frame(do.call(rbind, lapply(1:length(disparity), function(x) cbind(names(disparity)[x],t(disparity[[x]][[2]])))))
  names(plotmat)<-c("Name", "Disparity")
  plotmat$Age<-"Maastrichtian"
  plotmat$Age[grepl("Pal", plotmat$Name)]<-"Palaeocene"
  plotmat$Group<-sapply(plotmat$Name, function(x) strsplit(x, " ")[[1]][1])
  plotmat$Group[grepl("total", plotmat$Name)]<-"all"
  plotmat$Age<-factor(plotmat$Age, levels = unique(plotmat$Age))
  plotmat$Group<-factor(plotmat$Group, levels = unique(plotmat$Group))
  plotmat$Disparity<-as.numeric(plotmat$Disparity)
  
  ggplot(data=plotmat, aes(x=Group, y=Disparity, fill = Age))+
    geom_boxplot()+
    scale_fill_manual(values = c("#88C96F","#F9AA6F"), breaks = c("Maastrichtian","Palaeocene"), aesthetics = "fill") +
    xlab(label="Pollen groups")+
    ylab(label="Disparity (sum of variances)")+
    theme_bw()+
    theme(axis.ticks=element_blank())+
    labs(title = title)
}

###################

# Cycle through calculations with different data(sub)sets

results<-list("New all"=list(),
              "New South"=list(),
              "New North"=list(),
              "Lupia all"=list(),
              "Lupia KPg"=list())

for(dataset in names(results)){
  if(dataset %in% c("New all", "New South", "New North")){
    morphology <- read.csv("new morphology.csv", header = T, colClasses = "character")
    ranges<-read.csv("new ranges + groups.csv", header = T)

    # Remove c29 (uninformative) and other columns that aren't needed for this analysis
    morphology$c29<-NULL
    morphology$taxon_authorship<-NULL
    morphology$morphology_sources<-NULL
    morphology$notes<-NULL
    
    survivors<-switch(dataset,
                      "New all" = list("Maastrichtian" = which(ranges$Maastrichtian==1),
                                       "Boundary" = which(ranges$boundary==1),
                                       "Palaeocene" = which(ranges$Paleocene==1)),
                      "New South" = list("Maastrichtian" = which(ranges$pre_lat1==1),
                                         "Boundary" = c(),
                                         "Palaeocene" = which(ranges$post_lat1==1)),
                      "New North" = list("Maastrichtian" = which(ranges$pre_lat2==1),
                                         "Boundary" = c(),
                                         "Palaeocene" = which(ranges$post_lat2==1)))
  }else{
    morphology <- read.csv("Lupia1999 morphology.csv", header = T, colClasses = "character")
    ranges <- read.csv("Lupia1999 ranges + groups.csv", header = T)
    
    # remove characters ca1 - ca15 (last 7 columns)
    morphology <- morphology[,-c(32:38)]
    
    if(dataset=="Lupia KPg"){
      morphology <- morphology[ranges$Maastrichtian==1 | ranges$Paleocene==1,]
      ranges <- ranges[ranges$Maastrichtian==1 | ranges$Paleocene==1,]
    }
    
    survivors<-list("Maastrichtian" = which(ranges$Maastrichtian==1),
                    "Boundary" = c(),
                    "Palaeocene" = which(ranges$Paleocene==1))
  }
  
  # Set row names as taxon names and get rid of genus and species columns
  rownames(morphology) <- paste(morphology[,1],
                                morphology[,2],
                                sep = " ")
  
  rownames(ranges) <- paste(ranges[,1],
                            ranges[,2],
                            sep = " ")
  
  ranges<-ranges[,-c(1:2)]
  morphology<-morphology[,-c(1:2)]
  
  # Identify groups
  # Anemophilous: typically small, with minimal size variation, smooth-walled
  # or highly reduced surface ornamentation, thin-walled, with simple shape
  ranges$anemophilous<-as.integer(!(morphology$c7 == "2")
                                  & !(morphology$c9 == "2")
                                  & morphology$c16 %in% c("0", "1", "–", "?")
                                  & morphology$c21 %in% c("0", "–", "?")
                                  & morphology$c5=="0"
                                  & !morphology$c11 %in% c("2", "3", "4")
                                  & ranges$Triprojectacites==0)
  ranges$zoophilous<-as.integer(morphology$c7 == "2"
                                | morphology$c9 == "2"
                                | morphology$c16 %in% c("2", "3", "4", "5")
                                | !(morphology$c21 %in% c("0", "–", "?"))
                                | !morphology$c5=="0"
                                | morphology$c11 %in% c("2", "3", "4")
                                | ranges$Triprojectacites==1)
  ranges$ambivalent<-1-ranges$anemophilous-ranges$zoophilous
  
  groups<-list("All" = 1:nrow(ranges),
               "Colpate" = which(ranges$Colpate==1),
               "Triprojectate" = which(ranges$Triprojectacites==1),
               "Porate" = which(ranges$Porate==1),
               "Monosulcate" = which(ranges$Monosulcate==1),
               "Other" = which(ranges$Other==1),
               "Anemophilous" = which(ranges$anemophilous==1),
               "Zoophilous" = which(ranges$zoophilous==1),
               "Ambivalent" = which(ranges$ambivalent==1))
  
  grouplist<-as.list(outer(groups, survivors[-2], FUN=Vectorize(intersect)))[1:18]
  names(grouplist)<-outer(names(groups), names(survivors)[-2], paste)

  #############
  
  # Turn morphology dataframe into a cladistic matrix for Claddis
  # We need to replace the ? an – with NAs, and the characters A to E with 10 to 13
  # for this function to work
  morphology[morphology == "–"] <- NA
  morphology[morphology == "?"] <- NA
  morphology[morphology == "A"] <- 10
  morphology[morphology == "B"] <- 11
  morphology[morphology == "C"] <- 12
  morphology[morphology == "D"] <- 13
  morphology[morphology == "E"] <- 14
  
  lupia_cladmat <- build_cladistic_matrix(character_taxon_matrix = as.matrix(morphology),
                                          ordering = rep("unordered", ncol(morphology)))
  
  # Export matrix in nexus format, which means it can be read in more easily next time
  # (using Claddis's read_nexus_matrix function)
  write_nexus_matrix(cladistic_matrix = lupia_cladmat,
                     file_name = paste0("character matrix ", dataset, ".nex"))
  
  # Carry out ordination of the character matrix
  # This function is a wrapper, which both calculates a distance matrix and then
  # carries out a PCoA of that matrix
  lupia_pcoa <- ordinate_cladistic_matrix(cladistic_matrix = lupia_cladmat)
  
  # Extract PCoA scores and variance explained for plotting
  lupia_pcoa_scores <- lupia_pcoa$vectors
  lupia_pcoa_varexp <- lupia_pcoa$values$Rel_corr_eig*100 # convert to %
  
  #############
  
  # Disparity of K vs P taxa, using sums of variances of PCoA axis scores
  # Gives warnings for subsets with too few elements
  lupia_disp_sv<-disp_sv(lupia_pcoa_scores, grouplist, 1000)
  
  # Record run results
  results[[dataset]]<-list(cladmat=lupia_cladmat,
                           pcoa_scores=lupia_pcoa_scores,
                           pcoa_varexp=lupia_pcoa_varexp,
                           survivors=survivors,
                           groups=groups,
                           grouplist=grouplist,
                           disparity_sv=lupia_disp_sv)
  
  cat("'", dataset, "' batch done\n", sep = "")
}

# Save results
saveRDS(results, "disparity results.R")

# Summarise disparity results
dispdf<-do.call(rbind, lapply(1:length(results), \(x) summarise_disp(results[[x]]$disparity_sv, names(results)[x])))
write.csv(dispdf, "disparity summary.csv", row.names = FALSE)

# Calculate diversity/survival/extinction metrics
survival<-do.call(rbind, lapply(1:length(results), \(x) calc_survival(results[[x]]$grouplist, names(results)[x])))
write.csv(survival, "survival.csv", row.names = FALSE)

# Mann-Whitney U test on disparity changes
results_mannwu<-do.call(rbind, lapply(1:length(results), \(x) group_mannwu(results[[x]]$disparity_sv$disparity, names(results)[x])))
write.csv(results_mannwu, "mann-whitney u test.csv", row.names = FALSE)


############# Plotting figures

## Localities

occs<-read.csv("new occurrences KPg.csv", header = TRUE, fileEncoding = "UTF-8")

locs<-unique(occs[,4:5])
locs$bg<-ifelse(locs$latitude>48, "green3", "yellow")

chicxulub<-data.frame(lat=21.400000, long=-89.516667)
world<-map_data("world")

pdf("Figure 1 - Map.pdf", 5,5)
ggplot(world, aes(long, lat, group = group))+
  geom_polygon(fill="grey70", color="grey40", alpha=0.65, lwd=0.5)+
  coord_quickmap(xlim = c(-140, -75), ylim = c(17, 70))+
  geom_point(data = locs, aes(x=longitude, y=latitude, group = as.factor(latitude)), pch=21, color="black", bg=locs$bg, size=3)+
  geom_point(data = chicxulub, aes(x=long, y=lat, group = as.factor(long)), pch=21, color="black", bg="red", size=4)+
  geom_point(data = chicxulub, aes(x=long, y=lat, group = as.factor(long)), pch=21, color="black", bg="red", size=2)+
  labs(x="Longitude", y="Latitude")+
  theme(axis.ticks=element_blank())
dev.off()


## Taxonomic diversity

dispdf$Group<-sapply(dispdf$subsets, \(x) strsplit(x, " ")[[1]][[1]])
dispdf$Group<-factor(dispdf$Group, levels=names(groups), ordered = TRUE)
dispdf$Age<-sapply(dispdf$subsets, \(x) strsplit(x, " ")[[1]][[2]])

p1<-ggplot_div(dispdf, "New all", "A – New dataset")
p2<-ggplot_div(dispdf, "New South", "B – South bin")
p3<-ggplot_div(dispdf, "New North", "C – North bin")
p4<-ggplot_div(dispdf, "Lupia all", "D – Lupia (1999) dataset")

pdf("Figure 2 - Diversity.pdf", 10, 15)
grid.arrange(p1, p2, p3, p4, nrow=4)
dev.off()


## Extinction and survival

survival$Group<-sapply(survival$grouping, \(x) strsplit(x, " ")[[1]][[1]])
survival$Group<-factor(survival$Group, levels=names(groups), ordered = TRUE)
survival$Age<-sapply(survival$grouping, \(x) strsplit(x, " ")[[1]][[2]])

survival_survival<-survival[,names(survival) %in% c("batch", "Group", "Age", "survival_rate")]
survival_survival$Turnover<-"Survival"
names(survival_survival)[2]<-"value"
survival_extinction<-survival[,names(survival) %in% c("batch", "Group", "Age", "extinction_rate")]
survival_extinction$Turnover<-"Extinction"
names(survival_extinction)[2]<-"value"
survival_stack<-rbind(survival_survival, survival_extinction)
survival_stack<-survival_stack[!(survival_stack$Group %in% c("Other", "Ambivalent")),]

p1<-ggplot_survival(survival_stack, "New all", "A – New dataset")
p2<-ggplot_survival(survival_stack, "New South", "B – South bin")
p3<-ggplot_survival(survival_stack, "New North", "C – North bin")
p4<-ggplot_survival(survival_stack, "Lupia all", "D – Lupia (1999) dataset")

pdf("Figure 3 - Survivorship.pdf", 5, 15)
grid.arrange(p1, p2, p3, p4, nrow=4)
dev.off()


## PCoA

pdf(file = "Figure 4 - PCoAs.pdf", 10, 15)
layout(matrix(c(1:6), byrow = T, nrow = 3))
plot_pcoa_groups(results$`New all`$pcoa_scores, results$`New all`$pcoa_varexp, 1, 2, results$`New all`$groups, "A – New dataset pollen groups")
plot_survivors(results$`New all`$pcoa_scores, results$`New all`$pcoa_varexp, 1, 2, results$`New all`$survivors, "B – New dataset last occurrences")
plot_pcoa_groups(results$`Lupia all`$pcoa_scores, results$`Lupia all`$pcoa_varexp, 1, 2, results$`Lupia all`$groups, "C – Lupia (1999) whole dataset recalculated")
legend("topright", c("Monosulcates", "Tri-/polycolp(or)ates", "Triprojectacites", "Tri-/polyporates", "Others"), pch = c(15, 18, 18, 17, 19), bg="white", col = c("#A2FC3CFF", "grey70", "#455BCDFF", "#C42503FF", "grey90"), cex = 1.5, title = "Groups", box.col = "black")
plot_survivors(results$`Lupia all`$pcoa_scores, results$`Lupia all`$pcoa_varexp, 1, 2, results$`Lupia all`$survivors, "D – Lupia (1999) whole dataset last occurrences")
legend("topright", c("Maastrichtian", "Boundary clay", "Palaeocene"), pch = c(19, 19, 19), bg="white", col = c("grey20","green","#F9AA6F"), cex = 1.5, title = "Last occurrence", box.col = "black")
plot_pcoa_groups(results$`Lupia KPg`$pcoa_scores, results$`Lupia KPg`$pcoa_varexp, 1, 2, results$`Lupia KPg`$groups, "E – Lupia (1999) Maastrichtian–Palaeocene pollen groups")
plot_survivors(results$`Lupia KPg`$pcoa_scores, results$`Lupia KPg`$pcoa_varexp, 1, 2, results$`Lupia KPg`$survivors, "F – Lupia (1999) Maastrichtian–Palaeocene last occurrences")
dev.off()


## Disparity
p1<-ggplot_disp(results$`New all`$disparity_sv$disparity[!grepl("Ambivalent|Other", names(results$`New all`$disparity_sv$disparity))], "A – Total")
p2<-ggplot_disp(results$`New South`$disparity_sv$disparity[!grepl("Ambivalent|Other", names(results$`New South`$disparity_sv$disparity))], "B – South")
p3<-ggplot_disp(results$`New North`$disparity_sv$disparity[!grepl("Ambivalent|Other", names(results$`New North`$disparity_sv$disparity))], "C – North")

pdf(file = "Figure 5 - Disparities new.pdf", 8, 12)
grid.arrange(p1, p2, p3, nrow=3)
dev.off()


p1<-ggplot_disp(results$`Lupia all`$disparity_sv$disparity[!grepl("Ambivalent|Other", names(results$`Lupia all`$disparity_sv$disparity))], "A – Lupia total")
p2<-ggplot_disp(results$`Lupia KPg`$disparity_sv$disparity[!grepl("Ambivalent|Other", names(results$`Lupia KPg`$disparity_sv$disparity))], "B – Lupia Maastrichtian–Palaeocene")

pdf(file = "Figure 6 - Disparities Lupia.pdf", 8, 8)
grid.arrange(p1, p2, nrow=2)
dev.off()
