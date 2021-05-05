
library(adegenet); library(hierfstat)
library(parallel);	library(doParallel) #will load foreach
source("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/Fa_sample_funcs.R")
source("/home/user/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/Fa_sample_funcs.R")
colMax <- function(data) sapply(data, max, na.rm = TRUE)

#Get working directory
which_comp<-getwd()
	if (grepl("shoban",which_comp)) prefix_wd<-"C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/"
	if (grepl("shoban.DE",which_comp)) prefix_wd<-"C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/"
	if (grepl("home",which_comp)) prefix_wd<-"/home/user/Dropbox/Projects/IN_PROGRESS/"
setwd(prefix_wd)


#--------------------------------------------------------------------------------------------------
#	COMPARE EX SITU AND IN SITU NOW, AND THEN RESAMPLE WILD FOR MINIMUM SAMPLING/ OPTIMAL	#
#-------------------------------------------------------------------------------------------------

#####################################################
#													#
#	CALCULATE GEN DIV CAPTURE (% gen div in 		#
#	in all gardens and in each garden				#
#	AND MAKE .CSVs and R files for plots later		#
#													#
#####################################################

setwd(paste(prefix_wd,"Qhavardii_ex_situ/",sep=""))
gard_names<-read.csv("Key_to_Garden_POP.txt")[,2]

wild_files<-c("","_E","_W")
reg_names<-c("all","E","W")

#This will run over a loop of "Include all alleles (n_to_drop=0)" and "Include only alleles present in more than two copies (n_to_drop=2)"
for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-""
	if (n_to_drop==0) n_drop_file<-"_dr_0"

		#####################################################
		#													#
		#	Analysis by region- Total, East, West			#
		#	will output ex_vs_in_situ results csv by reg	#
		#													#
		#####################################################
		
	wild_results<-matrix(nrow=3,ncol=9+1)
	alleles_existing_by_sp<-matrix(nrow=3,ncol=9)
	#The loop below over 'reg' or 'regions' will go through three sets: analysis of all sample, only east samples and only west samples
	#The files are 'all wild pops plus all garden samples', 'all wild pops plus garden samples sourced from E', 'all wild pops plus garden samples sourced from W' 
	#The "wild_sets" will subset the wild reference to only the E and W for the reg=2 and reg=3 analysis
	#1-19 is east, 20-35 is west
	#The "garden_sets" idntify the populations that are garden samples
	#All three "by_garden" files have all wild pop'ns- will subset to east and west below
	#However the "_E" only has seedlings from East and "_W" only has seedlings from West
	#Note only four gardens have West seedlings
	
	wild_sets<-list(list(1,2:19,20:35),list(1,2:19),list(20,21:35))
	garden_sets<-(list(36:43,36:43,36:39))

	for (reg in 1:3){	
		Spp_tot_genind<-read.genepop(paste0("QH_total_garden_by_pop",wild_files[reg],".gen"),ncode=3)
		
		#This code compares the wild to various ex situ populations or all the ex situ merged
		#Just put in garden and wild population numbers
		wild_p<-unlist(wild_sets[[reg]]); garden_p<-unlist(garden_sets[[reg]])
		n_ind_W<-table(Spp_tot_genind@pop)[wild_p];  n_ind_G<-table(Spp_tot_genind@pop)[garden_p]; print(sum(n_ind_G))
		Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
		#Spp_tot_genind_sep<-seppop(Spp_tot_genind)
		alleles_cap<-colSums(Spp_tot_genpop[garden_p]@tab,na.rm=T)
		
		#Allele categories based only on wild populations (can look at all wild pop'ns or only one if you want)
		allele_cat_tot<-get.allele.cat(Spp_tot_genpop[wild_p], wild_sets[[reg]], 2, n_ind_W, glob_only=T,n_drop=n_to_drop)
		#This goes through each allele category and divides the number captured ex situ (alleles_cap) by the number of alleles existing (allele_cat_tot)
		for (i in 1:9) alleles_existing_by_sp[reg,i]<- (sum((allele_cat_tot[[i]])>0,na.rm=T))
		for (l in 1:length(allele_cat_tot)) wild_results[reg,l]<-round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
		
		#Basic genetic diversity stats alleles and Heterozygosity
		if(reg==1) print(paste("total allles", length(allele_cat_tot[[1]]), "alleles ex situ", sum(alleles_cap>0), "He E", mean(Hs(Spp_tot_genpop)[1:19]), "He W", mean(Hs(Spp_tot_genpop)[20:35]), "He gardens", mean(Hs(Spp_tot_genpop)[36:43])))
		
		wild_results[reg,10]<-sum(n_ind_G)
	}
	wild_results<-cbind(reg_names,wild_results)
		colnames(wild_results)<-c("region","global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare","num ind ex situ")
	write.csv(wild_results,file=paste("QH_ex_vs_in_situ_by_reg",n_drop_file,".csv",sep=""))
		
		#####################################################
		#													#
		#	Analysis by Garden- each of X gardens			#
		#	will output ex_vs_in_situ results csv by gard	#
		#													#
		#####################################################
		
	wild_results_by_g<-matrix(nrow=8,ncol=9+1)
	reg<-1
	#Go through each garden, thus garden_p (list of gardens) is replaced by gard i.e. just one garden
	for (gard in 36:43){	
		Spp_tot_genind<-read.genepop(paste0("QH_total_Garden_by_pop",wild_files[reg],".gen"),ncode=3)

		wild_p<-unlist(wild_sets[[reg]]); garden_p<-gard
		n_ind_W<-table(Spp_tot_genind@pop)[wild_p];  n_ind_G<-table(Spp_tot_genind@pop)[garden_p]; 
		Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
		alleles_cap<-colSums(Spp_tot_genpop[garden_p]@tab,na.rm=T)

		allele_cat_tot<-get.allele.cat(Spp_tot_genpop[wild_p], wild_sets[[reg]], 2, n_ind_W, glob_only=T,n_drop=n_to_drop)

		for (i in 1:9) alleles_existing_by_sp[reg,i]<- (sum((allele_cat_tot[[i]])>0,na.rm=T))
		
		for (l in 1:length(allele_cat_tot)) wild_results_by_g[gard-35,l]<-round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
		
		wild_results_by_g[gard-35,10]<-sum(n_ind_G)
	}

	wild_results_by_g<-cbind(as.character(gard_names),wild_results_by_g)
	colnames(wild_results_by_g)<-c("garden","global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare","num ind ex situ")
	write.csv(wild_results_by_g,file=paste("QH_ex_vs_in_situ_by_gard",n_drop_file,".csv",sep=""))
		
}	
	
		#####################################################
		#													#
		#	Optimal sampling by number trees				#
		#	Sample all possible sizes to create curve		#
		#	will output summ_results R file in subfolders	#
		#													#
		#####################################################
		
for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-""
	if (n_to_drop==0) n_drop_file<-"_dr_0"

	setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/")
	
	#for (reg in 1:3){ 
	reg<-1
		setwd(paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/Qhavardii_wild",wild_files[reg],sep=""))
		Spp_genind<-read.genepop(paste0("QH_wild",wild_files[reg],".gen"),ncode=3);	print(table(Spp_genind@pop))
		Spp_genpop<-genind2genpop(Spp_genind);	Spp_genind_sep<-seppop(Spp_genind)
			
		n_pops<-length(levels(Spp_genind@pop));		n_total_indivs<- length(Spp_genind@tab[,1])
		n_ind_p_pop<-table(Spp_genind@pop)
		allele_freqs<-colSums(Spp_genpop@tab)/(n_total_indivs*2)
		num_reps<-50000;		num_scen<-1
		alleles_existing_by_sp<-vector(length=9)
		
		list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
			
		allele_cat<-get.allele.cat(Spp_genpop, region_makeup=NULL, 2, n_ind_p_pop,n_drop=n_to_drop)	
		
		for (i in 1:9) alleles_existing_by_sp[i]<- (sum((allele_cat[[i]])>0,na.rm=T))
		#!!!RESULT!!!- HOW MANY ALLELES IN EACH CATEGORY	#proportion rare 
		alleles_existing_by_sp[5]/alleles_existing_by_sp[1]
		
		summ_results_tree<-array(dim=c(nrow(Spp_genind@tab)-1,11,num_scen,num_reps))
		for (nrep in 1:num_reps) {		#if not parallel comment IN
		#temp<-foreach(nrep=1:num_reps) %dopar% {	#if not parallel comment OUT
			alleles<-matrix(nrow=nrow(Spp_genind@tab)-1,ncol=length(allele_freqs))
			for (t in 2:(nrow(Spp_genind@tab)-1)){
				alleles<-colSums(Spp_genind@tab[sample(1:nrow(Spp_genind@tab), t),],na.rm=T)
				for (l in 1:length(allele_cat)) summ_results_tree[(t),(l+2),num_scen,nrep]<-sum(alleles[allele_cat[[l]]]>0, na.rm=T)
				}
		#		summ_results_tree[,,,nrep]	#if not parallel comment OUT
			}
		#summ_results_tree[,,1,]<-abind(temp,along=3)	#if not parallel comment OUT
		save(summ_results_tree,file=paste("summ_results_tree",n_drop_file,".R",sep=""))
	#}
}	
	
		#####################################################
		#													#
		#	Genetic distances- FST, clusters				#
		#													#
		#####################################################
	
	#FST
	reg<-1
	Spp_tot_genind<-read.genepop(paste0("QH_total_garden_for_FST",wild_files[reg],".gen"),ncode=3)
	m<-pairwise.fst(Spp_tot_genind)
	 makeSymm <- function(m) {
		m[upper.tri(m)] <- t(m)[upper.tri(m)]
		return(m)
	}
	a<-matrix(ncol=9,nrow=9)
	a[lower.tri(a)]<-m
	a<-makeSymm(a)
	rowMeans(a,na.rm=T);	mean(a[1,3:9],na.rm=T); mean(a[2,3:9],na.rm=T)	#FST East 0.0066, West 0.0159
	t.test(a[1,3:9],a[2,3:9],pair=T)	#p=0.049
		
	#Clusters
	dapc2 <- dapc(Spp_tot_genind,n.pca=50,n.da=50)
	pdf(file="garden_v_E_W_dapc.pdf")	
	scatter(dapc2,label=c("East","West",rep("G",8)),col=c("blue","red",rep("darkgrey",8)))
	dev.off()
	
	
	



#--------------------------------------------------------------------------------------------------
#	MAKE PLOTS AND DO STATISTICS	#
#-------------------------------------------------------------------------------------------------
		
library(adegenet);	library(hierfstat)

#Get working directory NOTE we are starting with IMLS Species then adding havardii so start in IMLS folder
which_comp<-getwd()
	if (grepl("shoban",which_comp)) prefix_wd<-"C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
	if (grepl("shoban.DE",which_comp)) prefix_wd<-"C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
	if (grepl("home",which_comp)) prefix_wd<-"/home/user/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
setwd(prefix_wd)
	
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare")


#####################################
#	COMPARE EX SITU OF HAVARDII TO OTHER SPECIES	
#	This figure is a plot of current conservation of genetic diversity in the current collection sizes
#	One point per species, and a log relationship across them, with havardii in focus
#	The first chunk of code creates plots of all allele types and over Reduced and Full Datasets
#	The second chunk of code is only the three alleles shown in figure in paper 
#####################################

##############################
#	First a set of plots for my own use, to see all results
##############################

for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-"";		if (n_to_drop==0) n_drop_file<-"_dr_0"
		
	wild_results_allsp<-read.csv(paste("ex_vs_in_situ",n_drop_file,".csv",sep=""))[,-1]
	wild_results<-data.frame(read.csv(file=paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/","QH_ex_vs_in_situ_by_reg",n_drop_file,".csv",sep=""))[,-1])
	names(wild_results_allsp)<-names(wild_results)
	wild_results_allsp<-rbind(data.frame(wild_results_allsp),data.frame(wild_results))
	
	pdf(file=paste("Qh_captured_ex_situ",n_drop_file,".pdf",sep=""),width=5,height=5)	
	for (i in 2:6){
	plot(wild_results_allsp[1:11,11],as.numeric(wild_results_allsp[1:11,i])*100,ylim=c(0,100), xlim=c(0,300),main="",xlab="number of plants in ex situ collections",ylab="percentage of alleles captured")
	text(as.numeric(wild_results_allsp[,11]),as.numeric(wild_results_allsp[,i])*100+5,substr(wild_results_allsp[,1],1,2),col=c(rep("darkgrey",11),"red","blue","black"),cex=c(rep(1,11),rep(1.5,3)))
	#Fitting data		
		d <- as.data.frame(list(plants=as.numeric(wild_results_allsp[1:11,11]), gendiv=as.numeric(wild_results_allsp[1:11,i])*100))
		mod <- lm(gendiv ~ I(log(plants)), d)		#Note: a square root relationship was tested and is ok but not as good as log
		new.df <- data.frame(plants=seq(10,250,by=1))
		out <- predict(mod, newdata = new.df)
		lines(unlist(c(new.df)), out, col = "grey", lwd=1.5)
		text(193,10,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black")
		points(wild_results_allsp[12:14,11],as.numeric(wild_results_allsp[12:14,i])*100,col=c("red","blue","black"),pch=19)
	}
	dev.off();		
} #close num to drop loop		


##############################
#	Now a plot for publication
##############################

pdf(file="Qh_captured_ex_situ_for_pub.pdf",width=9,height=4)	
	par(mfrow=c(1,3),mar=c(2,1.5,2,1),oma=c(4,5,2,1))
for (panel in 1:3){
	if (panel==1) {n_to_drop<-0; i<-2}
	if (panel==2) {n_to_drop<-2; i<-2}
	if (panel==3) {n_to_drop<-0; i<-5}
	if (n_to_drop==2) n_drop_file<-"";		if (n_to_drop==0) n_drop_file<-"_dr_0"
		
	wild_results_allsp<-read.csv(paste("ex_vs_in_situ",n_drop_file,".csv",sep=""))[,-1]
	wild_results<-data.frame(read.csv(file=paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/","QH_ex_vs_in_situ_by_reg",n_drop_file,".csv",sep=""))[,-1])
	names(wild_results_allsp)<-names(wild_results)
	wild_results_allsp<-rbind(data.frame(wild_results_allsp),data.frame(wild_results))
	
	plot(wild_results_allsp[1:11,11],as.numeric(wild_results_allsp[1:11,i])*100,ylim=c(0,100), xlim=c(0,300),main="",xlab="number of plants in ex situ collections",ylab="percentage of alleles captured",cex.axis=1.2)
	#show only the rare Quercus as letters- 7 to 9
	text(as.numeric(wild_results_allsp[c(7:9,12),11]),as.numeric(wild_results_allsp[c(7:9,12),i])*100+5,c(substr(wild_results_allsp[7:9,1],1,2),"Qh"),col=c(rep("darkgrey",3),"red","blue","black"),cex=c(rep(1.4,3),rep(1.6,3)))
	
	d <- as.data.frame(list(plants=as.numeric(wild_results_allsp[1:11,11]), gendiv=as.numeric(wild_results_allsp[1:11,i])*100))
	mod <- lm(gendiv ~ I(log(plants)), d)	
	new.df <- data.frame(plants=seq(10,250,by=1))
	out <- predict(mod, newdata = new.df)
	lines(unlist(c(new.df)), out, col = "grey", lwd=1.5)
	text(193,10,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black",cex=1.5)
	points(wild_results_allsp[12,11],as.numeric(wild_results_allsp[12,i])*100,col=c("red"),pch=19)
}	
mtext("Number of samples in ex situ collections currently",side=1,line=2,outer=T)
mtext("percent alleles captured",side=2,line=1.5,outer=T)
mtext("All alleles, Full Dataset                     All alleles, Reduced Dataset                    Low Frequency Alleles    ",side=3,line=0,outer=T)
dev.off()


#######################
#
#	Plot for OPTIMAL sampling- this plots the resampling output- every possible sampling size
#
####################### 
 
 
min_thresh<-95; y_lower<-90	# y_lower is the lower y axis limit for the plot
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare")

prefix_wd_imls<-"C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
folder_names<-c("Hhannerae","Hwaimeae","Masheii","Mpyramidata","Pekmanii","Psargentii","Qboyntonii","Qgeorgiana","Qoglethorpensis","Zdecumbens","Zlucayana","Qhavardii")	#	,"Qhavardii_E","Qhavardii_W") 

#calc_min will be a matrix to hold numbers for the minimum collection size thresholds for "sufficiency"
calc_min<-matrix(nrow=length(folder_names),ncol=length(list_allele_cat))
rownames(calc_min)<-c("H. w. hannerae","H. w. waimeae","M. ashei","M. pyramidata","P. ekmanii","P. sargentii","Q. boyntonii","Q. georgiana","Q. oglethorpensis","Z. decumbens","Z. lucayana", "Q. havardii")	#	, "Q. havardii East", "Q. havardii West")
colnames(calc_min)<-list_allele_cat


for (n_to_drop in c(2,0)){
	if (n_to_drop==0) n_drop_file<-"_dr_0" 
	if (n_to_drop==2) n_drop_file<-""
	pdf(file=paste("bytrees_overlay_merged_lowfreq_t_",n_drop_file,"2.pdf",sep=""),width=10, height=12)
	
par(mfrow=c(2,1))
	
for (i in c(3,6)){

if (i==3) { main_title<-"(A) All alleles"; x_upper_lim<-400 }
if (i==6) { main_title<-"(B) Low Frequency Alleles"; x_upper_lim<-150 }

	sp_colors<-c(rep("grey",11),"red")	#	,"purple","pink")	

	for (sp in 1:length(folder_names)){

		this_species<-folder_names[sp]
		if (sp<12) setwd(paste(prefix_wd_imls,this_species,sep=""))
		if (sp>=12) setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/Qhavardii_wild")

		load(file=paste("summ_results_tree",n_drop_file,".R",sep="")); num_reps<-length(summ_results_tree[1,1,1,])

		for (n in 1:num_reps) summ_results_tree[,,1,n]<-t(t(summ_results_tree[,,1,n])/summ_results_tree[length(summ_results_tree[,1,1,1]),,1,n])	
		#mean across reps
		all_mean<-apply(summ_results_tree[,,1,1:num_reps],c(1,2),mean,na.rm=T)*100
		#print(all_mean)
		
		if (sp==1) {
		plot(all_mean[,i],ylim=c(y_lower,100),type="l",lwd=4,xlim=c(0,x_upper_lim),xlab="number of plants (using simulated sampling)",
			ylab="percentage of genetic variation", cex.lab=1.2,col=sp_colors[sp],main=main_title)
		}	else {lines(all_mean[,i],lwd=4,col=sp_colors[sp])}
		abline(h=95,lty=2)
		legend(x_upper_lim*.7,93,c("11 rare species","Q. havardii"),col=c("grey","red"),lty=1,lwd=3)
		calc_min[sp,i-2]<-min(which(all_mean[,i]>min_thresh))
	}	#end species loop
}	
	write.csv(calc_min,file=paste("min_for_",min_thresh,n_drop_file,"2.csv",sep=""))	#make sure this works
	dev.off();		
#end to drop loop
}

########################
#Just to get number needed for 95% for havardii
##########################

setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/Qhavardii_wild")
load(file=paste("summ_results_tree",n_drop_file,".R",sep="")); num_reps<-1000

for (n in 1:num_reps) summ_results_tree[,,1,n]<-t(t(summ_results_tree[,,1,n])/summ_results_tree[length(summ_results_tree[,1,1,1]),,1,n])	
			
#mean across reps
all_mean<-apply(summ_results_tree[,i,1,1:num_reps],1,mean,na.rm=T)*100
print(all_mean)



#######################
#
#	LINEAR MODEL ON NUMBER PLANTS
#
####################### 

	setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Qhavardii_ex_situ/")

n_to_drop<-2; n_drop_file<-""
wild_results_by_g<-read.csv(file=paste("QH_ex_vs_in_situ_by_gard",n_drop_file,".csv",sep=""))[,-1]
colnames(wild_results_by_g)<-c("garden", "All", "Very common", "Common", "Low frequency", "Rare", "reg1", "loc1", "loc2", "loc 3", "num individuals")

pdf(file="linear_models_gardens.pdf",width=8,height=4)
par(mfrow=c(1,3),oma=c(4,4,2,2),mar=c(2,2,2,1))
for (i in c(2,5,6)){
	plot(wild_results_by_g[,11],as.numeric(wild_results_by_g[,i])*100,ylim=c(0,100), ylab="",xlab="",main=paste(colnames(wild_results_by_g)[i], "alleles"),cex.axis=1.3, pch=19, cex=1.4)
	#Fitting data		
		d <- as.data.frame(list(plants=wild_results_by_g[,11], gendiv=wild_results_by_g[,i]*100))
		#d<-rbind(c(0,0),d)
		if (i==2) mod <- lm(gendiv ~ I(sqrt(plants)), d)		
		if (i>2) mod <- lm(gendiv ~ I(log(plants)), d)		
		if (i==6) mod <- lm(gendiv ~ I((plants)), d)		
		new.df <- data.frame(plants=seq(5,200,by=1))
		out <- predict(mod, newdata = new.df)
		#plot(d,pch = 16, xlim = c(0,250), ylim = c(0,100))
		lines(unlist(c(new.df)), out, col = "red", lwd=1.5)
		text(100,10,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black",cex=1.3)
	}
	mtext("percent of alleles captured",line=2,side=2,outer=T)
	mtext("number of plants in garden",line=2,side=1,outer=T)
	dev.off()
	
	#plotting residuals
	par(mfrow=c(3,3))
for (i in c(2,5,6)){
	d <- as.data.frame(list(plants=wild_results_by_g[,11], gendiv=wild_results_by_g[,i]*100))		
	mod <- lm(gendiv ~ I(sqrt(plants)), d);	res <- residuals(mod)
	#qqnorm(res); qqline(res)
	plot(fitted(mod),res)
	mod <- lm(gendiv ~ I(log(plants)), d);	res <- residuals(mod)
	#qqnorm(res); qqline(res)
	plot(fitted(mod),res)
	mod <- lm(gendiv ~ I((plants)), d);	res <- residuals(mod)
	#qqnorm(res);  qqline(res)
	plot(fitted(mod),res)
	}
	
	
#######################
#
#	GETTTING NUM ACCESSIONS AND POPULATIONS
#
####################### 
	
qh_names<-read.csv("naming_temp.txt",sep="\t")
qh_names[qh_names[,5]=="The Morton Arboretum",3]
sort(unique(qh_names[qh_names[,5]=="Tulsa Botanic Garden",4]))
sort(unique(qh_names[qh_names[,5]=="Denver Botanical Garden",4]))
sort(unique(qh_names[qh_names[,5]=="Chicago Botanic Garden",4]))
sort(unique(qh_names[qh_names[,5]=="Boyce Thompson Arboretum (BTA)",3]))
sort(unique(qh_names[qh_names[,5]=="Lady Bird Johnson",3]))

	
#######################
#
#	LINEAR MODEL ON NUM PLANTS, ACCESSIONS, POPULATIONS SOURCED
#
####################### 


	all_all<-read.csv("all_temp_lm.txt",sep="\t")
pdf(file="linear_models_all_gardens.pdf",width=9,height=9)
par(mfrow=c(3,3),oma=c(4,4,2,3),mar=c(2,2,2,1))
for (i in c(2,5,6)){
for (j in c(7,9,8)){
	plot(all_all[,j],as.numeric(all_all[,i])*100,ylim=c(0,100), ylab="",xlab="",cex.axis=1.3, pch=19, cex=1.4)
	#Fitting data		
		d <- as.data.frame(list(plants=all_all[,j], gendiv=all_all[,i]*100))
		#d<-rbind(c(0,0),d)
		if (i==2) mod <- lm(gendiv ~ I(sqrt(plants)), d)		
		if (i>2) mod <- lm(gendiv ~ I(log(plants)), d)		
		if (i==6) mod <- lm(gendiv ~ I((plants)), d)		
		new.df <- data.frame(plants=seq(2,300,by=1))
		out <- predict(mod, newdata = new.df)
		#plot(d,pch = 16, xlim = c(0,250), ylim = c(0,100))
		lines(unlist(c(new.df)), out, col = "red", lwd=1.5)
		text(0.2*max(all_all[,j]),90,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black",cex=1.1)
	}
	mtext("percent of alleles captured",line=2,side=2,outer=T,cex=1.1)
	mtext("  number of plants in garden                 number of accessions               number of populations sourced",line=2,side=1,outer=T,cex=1.1)
	mtext("Globally Rare alleles                Globally Low frequency alleles                 All alleles                    ",line=1,side=4,outer=T,cex=1.1)

	}
	dev.off()
	

	
#######################
#
#	PLOT ALLELES CONTAINED IN EACH GARDEN
#
####################### 
	
	
pdf(file="by_garden_barplots.pdf",h=5,w=10)
	all_all<-read.csv("all_temp_lm.txt",sep="\t")
a<-barplot(t(as.matrix(all_all[,c(2,4:6)])),beside=T,names=c(paste("garden",LETTERS[1:8],"\nn=",all_all[1:8,7]),"all\ngardens"),col=c("grey21","grey42","grey66","grey93"),ylab="percent alleles conserved")
abline(h=.95,lty=2,col="grey36",lwd=1); abline(h=.7,lty=2,col="grey39")
legend(1,.90,c("all","common", "low frequency", "rare"),fill=c("grey21","grey42","grey63","grey93"),bg="white")
dev.off()
