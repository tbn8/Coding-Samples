#####Implement exact matching and greedy nearest neighbor matching on propensity scores, and test covariate balances.####

# load packages
library(tidyr)
library(plyr)
library(dplyr)
library(haven)
library(tableone)
library(reshape2)
library(ggplot2)

####### 1. MATCHING #######
# set globals
samplesize<-"full"
modelvars<-c("in_aco", "hrr", "M_75_84", "M_ge_85", "F_65_74", "F_75_84", "F_ge_85",
				"in_mdcd", "orig_disabled", "partd","has_ip_stay",
				"exceed_med_ip_days", "lt_p25", "bt_p75_p95", "ge_p95")
gtitle<-"Propensity Score Matching on All Variables + Exact Matching on Age-Gender-HRR"
modelname<-"exact_allvars_"

# read data
raw<-read_sas(paste0("../temp/unmatched_2018_",samplesize,".sas7bdat"))
keepvars<-c("BENE_ID", "in_aco", "hrr", "M_65_74", "M_75_84", "M_ge_85", "F_65_74", "F_75_84", "F_ge_85",
				"in_mdcd", "orig_disabled", "partd","has_ip_stay",
				"exceed_med_ip_days", "lt_p25", "bt_p75_p95", "ge_p95")
raw<-raw[,names(raw) %in% keepvars]
raw$hrr<-factor(raw$hrr) #code as categorical var

# propensity score calculation
Q<-glm(cbind(in_aco,1-in_aco) ~., family='binomial',data=raw[,names(raw) %in% modelvars])
ps<-predict(Q)
caliper<-(sd(ps*.25)/2)
raw<-transform(raw, ps=ps)
save(raw, file=paste0("../temp/raw_",modelname,samplesize,".Rdata"))

# split data by age-gender-HRR
agegrp<-c("M_65_74", "M_75_84", "M_ge_85", "F_65_74", "F_75_84", "F_ge_85")
raw$agesex<-apply(raw[,agegrp], 1, function(x) agegrp[which(x==1)])
demo<-split(raw, f=list(raw$hrr,raw$agesex), drop=TRUE)

# apply matching for each age-gender-HRR group
final<-lapply(demo, function(group) {

	# make treatment sample
	treatment<-group[group$in_aco==1,]
	N_treat<-dim(treatment)[1] #number of treatment

	if (N_treat==0) {
		NULL
	} else {
		# make control sample
		control<-group[group$in_aco==0,]
		numcol<-length(colnames(control)) #number of columns

		# matching
		msample<-data.frame(matrix(NA, nrow=N_treat, ncol=numcol))
		colnames(msample)<-colnames(control)
		set.seed(28)
		for (i in 1:N_treat) {
			temp<-control[control$ps < treatment$ps[i] + caliper & control$ps > treatment$ps[i] - caliper,]
			if (nrow(temp) > 0) {
				msample[i,]<-temp[(sample(1:nrow(temp), 1, replace=FALSE)),]
			}
		}
		# make crosswalk
		msample<-transform(msample,treatmentID=treatment$BENE_ID)
	}
})
# append all groups
final<-final[!sapply(final,is.null)]
msample<-do.call(rbind,final)
save(msample, file=paste0("../temp/matched_",modelname,samplesize,".Rdata"))

####### 2. TEST BALANCE #######
keepvars<-c("BENE_ID", "in_aco", "M_65_74", "M_75_84", "M_ge_85", "F_65_74", "F_75_84", "F_ge_85",
      "in_mdcd", "orig_disabled", "partd","has_ip_stay",
      "exceed_med_ip_days", "lt_p25", "bt_p75_p95", "ge_p95", "ps")

# create matched sample from crosswalk
control<-msample[,names(msample) %in% keepvars]
treatment<-raw[raw$in_aco==1,names(raw) %in% keepvars]
matched<-rbind(control,treatment)
head(matched)

# create random unmatched sample
set.seed(1234)
unmatched <- sample_n(raw, nrow(matched))

datasets <- c("matched", "unmatched")

# variables to test balance on
balvars <- c("M_65_74", "M_75_84", "M_ge_85", "F_65_74", "F_75_84", "F_ge_85",
      "in_mdcd", "orig_disabled", "partd","has_ip_stay",
      "exceed_med_ip_days", "lt_p25", "bt_p75_p95", "ge_p95", "ps")

newnames <- c("M 65-74", "M 75-84", "M >=85", "F 65-74", "F 75-84", "F >=85",
          "Medicaid", "Disabled", "Part D", "Has IP stay", "> Median IP stay",
          "Spend <25th", "Spend 75th-95th", "Spend >95th", "Propensity Score")

## construct standardized mean difference (SMD) tables
tabUnmatched <- CreateTableOne(vars=balvars, strata="in_aco", data = unmatched, test=FALSE)
tabMatched <- CreateTableOne(vars=balvars, strata="in_aco", data=matched, test=FALSE)

## construct a data frame containing variable names and SMD from all samples
dataPlot <- data.frame(variable  = rownames(ExtractSmd(tabUnmatched)),
                       Unmatched = as.numeric(ExtractSmd(tabUnmatched)),
                       Matched   = as.numeric(ExtractSmd(tabMatched)))

## create long-format data for ggplot2
dataPlotMelt <- melt(data          = dataPlot,
                     id.vars       = c("variable"),
                     variable.name = "Sample",
                     value.name    = "SMD")

## order factor levels in the same order
dataPlotMelt$variable <- factor(dataPlotMelt$variable,
                               levels = balvars)

## standardized mean difference plot
smd <- ggplot(data = dataPlotMelt,
       mapping = aes(x = variable, y = SMD, group = Sample, color = Sample, shape=Sample)) +
      geom_point(size=3,stroke=1) +
      scale_shape_manual(values=c("circle filled","triangle filled")) +
      geom_hline(yintercept = 0, color = "black", size = 0.5) +
      coord_flip() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_blank(),
            axis.text.y = element_text(size=12),
            axis.ticks.y = element_blank(),
            plot.title = element_text(size=16,face="bold"),
            legend.key = element_blank()) +
      labs(title = gtitle, x="", y="Absolute Standardized Mean Difference") +
      scale_x_discrete(labels=newnames)
ggsave(paste0("../output/smd_",modelname,samplesize,".png"),width=11,height=8.5)

## density plot for propensity scores
matched$status <- ifelse(matched$in_aco==1,"ACO","Non-ACO")
ps_density <- ggplot(data=matched) +
              stat_density(aes(x=ps,color=status),geom="line",position="identity") +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    plot.title = element_text(size=16,face="bold"),
                    legend.key = element_blank()) +
              labs(title=gtitle, x="Linear Probabilities",y="Density")
ggsave(paste0("../output/density_",modelname,samplesize,".png"),width=11,height=8.5)
