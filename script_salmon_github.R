library(brms)
library(ggridges)
library(tidyverse)
library(scales)
library(viridis)
library(ggstance)
library(rstan)
#######IMPORTANT#####
#####You will need to download the C++ compiler outside of R first###
#####Nothing will work if it is not installed first###
#####Follow these instructions to download the C++ compiler and install RStan (which runs the Hamiltonian Monte Carlo)
#####https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows



####model#####
#priors for intercept chosen to cover mean (39, log=3.6) with a standard deviation that would cover the max at 2sd away from log(3.6) flux of 
#emergence measured from g_m_yr studies in Gratton and Vander Zanden 2009. File=gratton_data.csv.
data$gcmyr01<-0.01+data$gcmyr #add 0.01 to flux

m2<-brm(gcmyr01~(1|River)+(1|site/year),data=data,family=Gamma(link="log"),
        prior=c(prior(normal(3.6,3),class="Intercept"),
                prior(cauchy(0,1),class="sd")),
        chains=4,iter=5000)

m2 #summarize model
pp_check(m2,type="boxplot") #posterior predictive checks
m2post<-posterior_samples(m2) #sample joint posterior

###Prediction######
newdata<-data.frame(River="new", #create new data frame to predict for. 
                    site="new",
                    year=2019)

m2pred<-data.frame(predict(m2,type="response",newdata=newdata,re_formula=~(1|River)+(1|site/year),
                           allow_new_levels  = TRUE,summary=FALSE)) #simulate predictions from the posterior predictive distribution
m2preda<-gather(m2pred) #change format of predictions from wide to long (to merge later with fitted predictions for plotting)
m2preda$key<-"predicted" #id for predicted estimates
m2preda$color<-"predicted" #id for predicted estimates

###Fitted estimates ####
#estimates of mean flux across sample sites
m2int<-data.frame(exp(m2post$b_Intercept)) #fitted distribution of mean flux (intercept) across sites
m2int<-gather(m2int) #change format from wide to long (to merge with predictions)
m2int$key<-"fitted" #id for fitted estimates
m2int$color<-"fitted" #id for fitted estimates


###Convert flux (g/m/yr) to deposition (g/m2/yr) estimates######
m2plot2<-rbind(m2int,m2preda) #combine fitted and predicted estimates

m2plot2$prop_males<-0.55 #from supplementary of Walters et al. (2018)
m2plot2$exuvia_f_gcmyr<-m2plot2$value*0.45*.13 #female exuviae. 45% of population are female. 13% of female weight is exuvia
m2plot2$exuvia_m_gcmyr<-m2plot2$value*0.55*0.21 #male exuviae. 55% of populations are male. 21% of male weight is exuvia
m2plot2$exuvia_tot<-m2plot2$exuvia_f_gcmyr+m2plot2$exuvia_m_gcmyr #sum of female and male exuviae fractions. Add this to adult estiamtes

m2plot2$adult_C<-m2plot2$value-m2plot2$exuvia_tot #adult carbon deposition
m2plot2$adult_P<-m2plot2$adult_C/124 #adult P deposition. assumes a C:P of 124 from Elser et al. 2000 for aquatic insect herbivores
m2plot2$adult_N<-m2plot2$adult_C/6.3 #adult N deposition. assumes a C:P of 6.3 from Elser et al. 2000 for aquatic insect herbivores
m2plot2$exuvia_C<-m2plot2$exuvia_tot*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018
m2plot2$exuvia_P<-m2plot2$exuvia_C/124*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018
m2plot2$exuvia_N<-m2plot2$exuvia_C/6.3*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018

#CARBON deposition
m2plot2$flux_50c<-m2plot2$adult_C/2/2.7 #divide flux by 2 to get 50%, then by 2.7 to get estiamte of 100% deposition.
m2plot2$dep_100c<-m2plot2$flux_50c+m2plot2$exuvia_C #assume that 100% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_25c<-m2plot2$flux_50c*0.25+m2plot2$exuvia_C #assume that 25% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_50c<-m2plot2$flux_50c*0.50+m2plot2$exuvia_C #assume that 50% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_75c<-m2plot2$flux_50c*0.75+m2plot2$exuvia_C #assume that 75% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$flux_50c<-NULL #remove column - no longer needed.

#NITROGEN deposition
m2plot2$flux_50n<-m2plot2$adult_N/2/2.7 #divide flux by 2 to get 50%, then by 2.7 to get estiamte of 100% deposition.
m2plot2$dep_100N<-m2plot2$flux_50n+m2plot2$exuvia_N #assume that 100% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_25N<-m2plot2$flux_50n*0.25+m2plot2$exuvia_N #assume that 25% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_50N<-m2plot2$flux_50n*0.50+m2plot2$exuvia_N #assume that 50% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_75N<-m2plot2$flux_50n*0.75+m2plot2$exuvia_N #assume that 75% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$flux_50n<-NULL

#PHOSPHOROUS deposition
m2plot2$flux_50p<-m2plot2$adult_P/2/2.7 #divide flux by 2 to get 50%, then by 2.7 to get estiamte of 100% deposition.
m2plot2$dep_100P<-m2plot2$flux_50p+m2plot2$exuvia_P #assume that 100% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_25P<-m2plot2$flux_50p*0.25+m2plot2$exuvia_P #assume that 25% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_50P<-m2plot2$flux_50p*0.50+m2plot2$exuvia_P #assume that 50% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$dep_75P<-m2plot2$flux_50p*0.75+m2plot2$exuvia_P #assume that 75% of available salmonflies (50% of adults + 100% of exuviae) are deposited to detrital pools
m2plot2$flux_50p<-NULL



##
#####FIGURES##########
##

#####Figure 2 - raw data plotted against benchmarks######
#add predictions of flux and deposition to m2pred
data$exuvia_f_gcmyr<-data$gcmyr*(1-data$sex)*.13 #female exuviae. 45% of population are female. 13% of female weight is exuvia
data$exuvia_m_gcmyr<-data$gcmyr*data$sex*0.21 #male exurviae. 55% of populations are male. 21% of male weight is exuvia
data$exuvia_tot<-data$exuvia_f_gcmyr+data$exuvia_m_gcmyr

data$adult_C<-data$gcmyr-data$exuvia_tot
data$adult_P<-data$adult_C/124 #assumes a C:P of 124 from Elser et al. 2000 for aquatic insect herbivores
data$adult_N<-data$adult_C/6.3 #assumes a C:P of 6.3 from Elser et al. 2000 for aquatic insect herbivores
data$exuvia_C<-data$exuvia_tot*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018
data$exuvia_P<-data$exuvia_C/124*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018
data$exuvia_N<-data$exuvia_C/6.3*0.65 #assumes that exuviae have 65% of the C, N, and P of adults via Callaham and Whiles 2001, and also via email communication with Callaham on 6/14/2018

data$flux_50c<-data$adult_C/2/2.7
data$dep_100c<-data$flux_50c+data$exuvia_C #amount of flux within 2.7 m of shore on a per m2 basis, assuming 50% falls within 2.7 m
data$dep_25c<-data$flux_50c*0.25+data$exuvia_C #assume that 25% of flux is deposited to detrital pools
data$dep_50c<-data$flux_50c*0.50+data$exuvia_C #assume that 50% of flux is deposited to detrital pools
data$dep_75c<-data$flux_50c*0.75+data$exuvia_C #assume that 75% of flux is deposited to detrital pools
data$flux_50c<-NULL

data$flux_50n<-data$adult_N/2/2.7
data$dep_100N<-data$flux_50n+data$exuvia_N #amount of flux within 2.7 m of shore on a per m2 basis, assuming 50% falls within 2.7 m
data$dep_25N<-data$flux_50n*0.25+data$exuvia_N #assume that 25% of flux is deposited to detrital pools
data$dep_50N<-data$flux_50n*0.50+data$exuvia_N #assume that 50% of flux is deposited to detrital pools
data$dep_75N<-data$flux_50n*0.75+data$exuvia_N #assume that 75% of flux is deposited to detrital pools
data$flux_50n<-NULL

data$flux_50p<-data$adult_P/2/2.7
data$dep_100P<-data$flux_50p+data$exuvia_P #amount of flux within 2.7 m of shore on a per m2 basis, assuming 50% falls within 2.7 m
data$dep_25P<-data$flux_50p*0.25+data$exuvia_P #assume that 25% of flux is deposited to detrital pools
data$dep_50P<-data$flux_50p*0.50+data$exuvia_P #assume that 50% of flux is deposited to detrital pools
data$dep_75P<-data$flux_50p*0.75+data$exuvia_P #assume that 75% of flux is deposited to detrital pools
data$flux_50p<-NULL

###format data for plotting
datag<-gather(data,est,gm2yr,"dep_100c":"dep_75P") #create new column with different deposition estimates (from 25-100%)
datag$element<-ifelse(grepl("c",datag$est),"carbon",
                      ifelse(grepl("P",datag$est),"phosphorous","nitrogen"))
datag$amount<-ifelse(grepl("100",datag$est),"100%",
                     ifelse(grepl("75",datag$est),"75%",
                            ifelse(grepl("50",datag$est),"50%","25%")))

datag$amount<-factor(datag$amount,levels=c("100%","75%","50%","25%"))

datag$atm_dep<-as.numeric(ifelse(datag$element=="nitrogen","0.35",""))
datag$atm_dep_upper<-as.numeric(ifelse(datag$element=="nitrogen","0.75",""))
datag$atm_dep_lower<-as.numeric(ifelse(datag$element=="nitrogen","0.1",""))



figure_2<-ggplot()+
  geom_point(data=datag,aes(x=reorder(site2,-gm2yr),y=gm2yr,color=element),
             size=1.5,alpha=0.9)+  
  geom_hline(data=datag,aes(yintercept=atm_dep),color="gray50",size=1.2)+
  geom_hline(data=Burns_deposition_data,aes(yintercept=tot_N_gm2yr),alpha=0.15)+
  geom_hline(data=terr_prod_Bartrons,aes(yintercept=terr_gm2yr),alpha=0.3)+
  geom_hline(data=P_dep_Tipping,aes(yintercept=terrP_gm2yr),alpha=0.3)+
  #geom_segment(data=datag,aes(x=reorder(xline,-gm2yr),xend=reorder(xend_line,-gm2yr),y=yline,yend=yend_line))+
  facet_grid(element~amount,scales="free")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        text=element_text(size=12),
        axis.title.y=element_text(size=12))+
  guides(color=FALSE)+
  #scale_y_log10()+
  scale_color_viridis(discrete=TRUE,option="E",end=0.3)+
  ylab(expression(paste("Potential detrital deposition of salmonflies (g/m"^2,"/y)")))+
  xlab("Site")

figure_2

ggsave(figure_2,file="figure_2.jpg",dpi=500,width=6,height=4.5,units="in")


######Figure 3#####
####format estimates for plotting
m2plot2<-gather(m2plot2,category,gm2yr,"dep_100c":"dep_75P")
m2plot2$amount<-ifelse(grepl("100",m2plot2$category),"100%",
                        ifelse(grepl("75",m2plot2$category),"75%",
                               ifelse(grepl("50",m2plot2$category),"50%","25%")))
m2plot2$amount<-factor(m2plot2$amount,levels=c("100%","75%","50%","25%"))
m2plot2$key<-factor(m2plot2$key,levels=c("predicted","fitted"))
m2plot2$element<-ifelse(grepl("c",m2plot2$category),"carbon",
                        ifelse(grepl("N",m2plot2$category),"nitrogen","phosphorous"))

figure_3<-ggplot(data=m2plot2,aes(x=gm2yr,y=key,fill=element))+
  geom_density_ridges2(alpha=0.9)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1,10),labels=c("0.0001","0.001","0.01","0.1","1","10"))+
  facet_grid(amount~.)+
  scale_fill_manual(values=c("#2e75b6","grey60","grey80"))+
  coord_cartesian(xlim=c(0.00002,10))+
  annotation_logticks(base=10,sides="b")+
  theme_classic()+
  xlab(expression(paste("Deposotion of salmonflies (g/m"^2,"/y)")))+
  theme(axis.title.y=element_blank())

figure_3
ggsave(figure_3,file="figure_3.jpg",dpi=500,width=5.5,height=7)  

####summary_stats and tables#####
#####data for Table 1####
summary_table2<-m2plot2%>%
  group_by(element,amount,key)%>%
  summarise('low95'=quantile(gm2yr,probs=0.025),
            'med'=quantile(gm2yr,probs=0.5),
            'high95'=quantile(gm2yr,probs=0.975),
            'mean'=mean(gm2yr),
            'sd'=sd(gm2yr))
write.csv(summary_table2,file="summary_table2.csv")




#######Supplementary Material#####


#model with wider prior for sd of intercept 
m2b<-brm(gcmyr01~(1|River)+(1|site/year),data=data,family=Gamma(link="log"),
         prior=c(prior(normal(3.6,6),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=5000)

m2b #did not converge

#model with smaller mean for indercept
m2c<-brm(gcmyr01~(1|River)+(1|site/year),data=data,family=Gamma(link="log"),
         prior=c(prior(normal(1,3),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=5000)

m2c

#extract posteriors
m2_sum<-data.frame(posterior_summary(m2))
m2b_sum<-data.frame(posterior_summary(m2b))
m2c_sum<-data.frame(posterior_summary(m2c))

#id each model
m2_sum$alternate_priors<-"N(3.6,3) model used"
m2b_sum$alternate_priors<-"N(3.6,6) sd X 2 (did not converge)"
m2c_sum$alternate_priors<-"N(1,3) reduce mean"

#combine estimates from the three models
m2_compare<-rbind(m2_sum,m2b_sum,m2c_sum)
m2_compare$parameter<-rownames(m2_compare)
m2_compare$effect<-ifelse(grepl("r_",m2_compare$parameter),"re_offset",
                          ifelse(grepl("sd",m2_compare$parameter),"group level_sd",
                                 ifelse(grepl("nterc",m2_compare$parameter),"Intercept",
                                        ifelse(grepl("shape",m2_compare$parameter),"family specific","lp"))))
m2_compare$alternate_priors<-factor(m2_compare$alternate_priors,levels=c("N(3.6,3) model used",
                                                                         "N(1,3) reduce mean",
                                                                         "N(3.6,6) sd X 2 (did not converge)"))


####Figure S1 - Prior vs posterior plots######

prior_int<-data.frame(rnorm(10000,3.6,3))  #intercept prior
post_int<-data.frame(rnorm(10000,0.9,0.42))  #intercept posterior (will vary with each run, but should be close to these values)

prior_int<-gather(prior_int)
post_int<-gather(post_int)

post_int$key<-"posterior intercept"
prior_int$key<-"prior intercept"
post_int$color<-"posterior"
prior_int$color<-"prior"

post_int$facet<-"intercept"
prior_int$facet<-"intercept"

post_priorm2<-rbind(post_int,prior_int)

figure_s1<-ggplot(post_priorm2,aes(x=value,fill=color))+
  geom_density(alpha=.5)+
  #facet_wrap(~facet,scales="free")+
  #coord_cartesian(xlim=c(0,50))+
  scale_fill_manual(values=c("grey50","blue"))+
  theme_classic()+
  ylab("Probability density")+
  theme(legend.title=element_blank())+
  xlab(expression(paste("Intercept (log salmonfly flux g/m/y)")))

figure_s1

ggsave(figure_s1,file="figure_s1.jpg",dpi=500,width=7,height=6,units="in")

#####Figure S2######

library(gridExtra)

p1<-ggplot(subset(m2_compare,effect!="re_offset"&effect!="lp"),aes(x=alternate_priors,y=Estimate,ymax=Q97.5,ymin=Q2.5,color=alternate_priors))+
  geom_pointrange(position=position_dodge(width=0.5))+
  facet_wrap(~effect,scales="free")+
  theme_classic()+
  scale_color_grey()+
  guides(color=FALSE)+
  theme(axis.text.x=element_text(angle=45,hjust=1),,
        text=element_text(size=15))

p2<-ggplot(subset(m2_compare,effect=="re_offset"),aes(x=parameter,y=Estimate,ymax=Q97.5,ymin=Q2.5,color=alternate_priors,fill=alternate_priors,group=parameter))+
  geom_pointrange(position=position_dodge(width=0.01),size=0.2)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        text=element_text(size=15))+
  scale_color_grey()+
  xlab("random effect offsets for each site, river, and year (labels removed for clarity)")

figure_s2<-grid.arrange(p2,p1,ncol=1)
ggsave(figure_s2,file="figure_s2.jpg",width=11,height=9)




