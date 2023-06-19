library(here)
library(ggplot2); 
library(cmdstanr)
library(samEst)
library(cowplot)
library(xtable)

mytheme = list(
  theme_classic(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

stock_info<- read.csv(here::here("data/psalmon_info.csv"))
stock_dat<-read.csv(here::here("data/psalmon_sr.csv"))


#eg. for Cowichan chinook

srhar<-subset(stock_dat,stock=="Harrison"&species=="Chinook")

hardat<-data.frame(stock=paste(srhar$stock,srhar$species),
                   by=srhar$broodyear,
                   S=srhar$spawners,
                   R=srhar$recruits,
                   logRS=log(srhar$recruits/srhar$spawners))

dirpr<-matrix(c(2,1,1,2),2,2)
psimple <- ricker_TMB(data=hardat) #note you might get an error here - other models might still fit though
pac<-ricker_TMB(data=hardat, AC=TRUE)
ptva<- ricker_rw_TMB(data=hardat,tv.par="a")
ptvb <- ricker_rw_TMB(data=hardat,tv.par="b",sig_p_sd=1)
ptvab <- ricker_rw_TMB(data=hardat,tv.par="both",sig_p_sd=1)
phmma <- ricker_hmm_TMB(data=hardat, tv.par='a',dirichlet_prior=dirpr)
phmmb <- ricker_hmm_TMB(data=hardat, tv.par='b',dirichlet_prior=dirpr)
phmm <- ricker_hmm_TMB(data=hardat, tv.par='both',dirichlet_prior=dirpr)


cphar<-c(psimple$conv_problem, 
pac$conv_problem, 
ptva$conv_problem, 
ptvb$conv_problem, 
ptvab$conv_problem, 
phmma$conv_problem, 
phmmb$conv_problem, 
phmm$conv_problem) 

#-------------------------------------------------
#predictions
x_new=seq(0,max(hardat$S)*1.5,length.out=200)
by_q=round(quantile(hardat$by,seq(0,1,by=0.1)))
#static
pred_simple=data.frame(pred=exp(psimple$alpha-psimple$beta*x_new)*x_new,
                       x_new=x_new,
                       gr=rep(.1,length(x_new)),
                       model="simple")
pred_ar1=data.frame(pred=exp(pac$alpha-pac$beta*x_new)*x_new,
                    x_new=x_new,
                    gr=rep(.1,length(x_new)),
                    model="ar1")

rsrwa=list()
for(n in 1:length(by_q)){
  rsrwa[[n]]=exp(ptva$alpha[match(by_q[n],hardat$by)]-ptva$beta*x_new)*x_new
}
pred_rwa=data.frame(pred=do.call(c,rsrwa),
                    x_new=x_new,
                    gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                    model="rwa")
alpha_ptva=data.frame(by=hardat$by,value=ptva$alpha,param="alpha")

rsrwb=list()
for(n in 1:length(by_q)){
  rsrwb[[n]]=exp(ptvb$alpha-ptvb$beta[match(by_q[n],hardat$by)]*x_new)*x_new
}
beta_ptvb=data.frame(by=hardat$by,value=ptvb$beta,param="beta")
pred_rwb=data.frame(pred=do.call(c,rsrwb),
                    x_new=x_new, 
                    gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                    model="rwb")


rsrwab=list()
for(n in 1:length(by_q)){
  rsrwab[[n]]=exp(ptvab$alpha[match(by_q[n],hardat$by)]-ptvab$beta[match(by_q[n],hardat$by)]*x_new)*x_new
}
beta_ptvab=data.frame(by=hardat$by,value=c(ptvab$alpha,ptvab$beta),param=rep(c("alpha","beta"),each=length(hardat$by)))

pred_rwab=data.frame(pred=do.call(c,rsrwab),
                     x_new=x_new,
                     gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                     model="rwab"
)


rshmma1=exp(phmma$alpha[1]-phmma$beta*x_new)*x_new
rshmma2=exp(phmma$alpha[2]-phmma$beta*x_new)*x_new
pred_hmma=data.frame(pred=c(rshmma1,rshmma2),
                     x_new=x_new,
                     gr=rep(c(.1,1),each=length(x_new)),
                     model="HMMa")

gamma_df_hmma=data.frame(by=hardat$by,gamma=phmma$probregime[2,])

rshmmb1=exp(phmmb$alpha-phmmb$beta[1]*x_new)*x_new
rshmmb2=exp(phmmb$alpha-phmmb$beta[2]*x_new)*x_new
pred_hmmb=data.frame(pred=c(rshmmb1,rshmmb2),
                     x_new=x_new,
                     gr=rep(c(1,.1),each=length(x_new)),
                     model="HMMb")

gamma_df_hmmb=data.frame(by=hardat$by,gamma=phmmb$probregime[2,])



rshmmab1=exp(phmm$alpha[1]-phmm$beta[1]*x_new)*x_new
rshmmab2=exp(phmm$alpha[2]-phmm$beta[2]*x_new)*x_new
pred_hmmab=data.frame(pred=c(rshmmab1,rshmmab2),
                      x_new=x_new,
                      gr=rep(c(.1,1),each=length(x_new)),
                      model="HMMab")

gamma_df_hmmab=data.frame(by=hardat$by,gamma=phmm$probregime[2,])

pred_har<-rbind(pred_simple,pred_ar1,pred_rwa,pred_rwb,pred_rwab,pred_hmma,pred_hmmb,pred_hmmab)

#--------------------------------------------
#prediction plot
pred_har$model<-factor(pred_har$model, levels=c("rwa", "rwb", "rwab","HMMa","HMMb", "HMMab","simple", "ar1") )

p<-ggplot(pred_har)+
  geom_line(aes(x=x_new,y=pred,colour = as.factor(gr)),linewidth=1.5)+
  geom_point(data=hardat, aes(x=S, y=R,alpha=by),size=1.5) +
  scale_colour_viridis_d(name='time')+
  facet_wrap(~model)+
  xlab("Spawners") + 
  ylab("Recruits")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
p
ggsave("tex/predhar.png",plot=p,width = 10, height = 7)

#==========================================
#parameter trends

dfpar<-data.frame(by=rep(hardat$by,16),
                  value=c(
                    rep(psimple$alpha,length(hardat$by)),
                    rep(pac$alpha,length(hardat$by)),
                    ptva$alpha,
                    rep(ptvb$alpha,length(hardat$by)),
                    ptvab$alpha,
                    phmma$alpha[phmma$regime],
                    rep(phmmb$alpha,length(hardat$by)),
                    phmm$alpha[phmm$regime],
                    1/rep(psimple$beta,length(hardat$by)),
                    1/rep(pac$beta,length(hardat$by)),
                          1/rep(ptva$beta,length(hardat$by)),
                          1/ptvb$beta,
                          1/ptvab$beta,
                          1/rep(phmma$beta,length(hardat$by)),
                          1/phmmb$beta[phmmb$regime],
                          1/phmm$beta[phmm$regime]),
                    param=rep(c("log_a","Smax"), each=length(hardat$by)*8),
                    model=factor(rep(rep(c("simple","AR1","RWa","RWb","RWab",
                                "Shifta","Shiftb","Shiftab"),2),each=length(hardat$by)),
                                levels=c("simple","AR1","RWa","RWb","RWab",
                                         "Shifta","Shiftb","Shiftab")))




pp<-ggplot(dfpar) +   
  geom_line(aes(x=by,y=value,color=by ),
            linewidth=1.4) +  
  mytheme+ 
  ylab("parameters") +
  xlab("year") +
  facet_grid(param~model, scales="free_y")+
  scale_color_viridis_c(begin=.1, end=.8, guide="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
pp
ggsave("tex/paramhar.png",plot=pp,width = 12, height = 5)



#model selection criteria
#lfo
#lfo
    lfostatic <- tmb_mod_lfo_cv(data=hardat,model='static', L=15)
    lfoac <- tmb_mod_lfo_cv(data=hardat,model='staticAC', L=15)
    lfoalpha <- tmb_mod_lfo_cv(data=hardat,model='rw_a', siglfo="obs", L=15)
    lfobeta <- tmb_mod_lfo_cv(data=hardat,model='rw_b', siglfo="obs", L=15)
    lfoalphabeta <- tmb_mod_lfo_cv(data=hardat,model='rw_both', siglfo="obs", L=15)
    lfohmma <- tmb_mod_lfo_cv(data=hardat,model='HMM_a', L=15)
    lfohmmb <- tmb_mod_lfo_cv(data=hardat,model='HMM_b', L=15)
    lfohmm <- tmb_mod_lfo_cv(data=hardat,model='HMM', L=15)
 lfores<-c(sum(lfostatic$lastparam), 
        sum(lfoac$lastparam), 
        sum(lfoalpha$lastparam), 
        sum(lfobeta$lastparam), 
        sum(lfoalphabeta$lastparam), 
        sum(lfohmma$lastregime_pick),
        sum(lfohmmb$lastregime_pick),
        sum(lfohmm$lastregime_pick))
                           


df<-data.frame(model=c("simple","AR1","RW a","RW b","RW a and b",
                    "Shift a","Shift b","Shift a and b"),
           AIC=round(c(psimple$AICc, pac$AICc,
                 ptva$AICc,ptvb$AICc, NA,
                 phmma$AICc,phmmb$AICc, phmm$AICc)-min(
                c(psimple$AICc, pac$AICc,
                     ptva$AICc,ptvb$AICc, ptvab$AICc,
                 phmma$AICc,phmmb$AI, ptvab$AICc,
                     phmma$AICc,phmmb$AICc, phmm$AICc)),2),
           BIC=round(c(psimple$BIC, pac$BIC,
                 ptva$BIC,ptvb$BIC, NA,
                 phmma$BIC,phmmb$BIC, phmm$BIC)-
             min(c(psimple$BIC, pac$BIC,
                   ptva$BIC,ptvb$BIC, ptvab$BIC,
                   phmma$BIC,phmmb$BIC, phmm$BIC)),2),
             LFO=round(c(sum(lfostatic$lastparam), 
                       sum(lfoac$lastparam), 
                       sum(lfoalpha$lastparam), 
                       sum(lfobeta$lastparam), 
                       sum(lfoalphabeta$lastparam), 
                       sum(lfohmma$lastregime_pick),
                       sum(lfohmmb$lastregime_pick),
                       sum(lfohmm$lastregime_pick)),2))
             

df$AIC = ifelse(df$AIC ==0, paste0("\\textbf{", df$AIC, "}"), df$AIC)
df$BIC = ifelse(df$BIC ==0, paste0("\\textbf{", df$BIC, "}"), df$BIC)

df$LFO = ifelse(df$LFO ==max(df$LFO, na.rm=T), paste0("\\textbf{", df$LFO, "}"), df$LFO)



print(xtable(df,caption="Parameter trajectories for stationary 
  and time varying models considered for Harrison Chinook. 
  Minimum AIC and BIC value indicate preferred model, while 
  maximum LFO indicate preferred model. Missing values 
  indicate lack of model convergence",
label="harcrit"), caption.placement="top", file="tex/harcrit.tex",
sanitize.text.function = identity,NA.string = "-")



#===================================
#cow
dirpr<-matrix(c(2,1,1,2),2,2)
srcow<-subset(stock_dat,stock=="Cowichan"&species=="Chinook")


cowdat<-data.frame(stock=paste(srcow$stock,srcow$species),
                   by=srcow$broodyear,
                   S=srcow$spawners,
                   R=srcow$recruits,
                   logRS=log(srcow$recruits/srcow$spawners))


gsimple <- ricker_TMB(data=cowdat) 
gac<-ricker_TMB(data=cowdat, AC=TRUE)
gtva<- ricker_rw_TMB(data=cowdat,tv.par="a")
gtvb <- ricker_rw_TMB(data=cowdat,tv.par="b",sig_p_sd=1)
gtvab <- ricker_rw_TMB(data=cowdat,tv.par="both",sig_p_sd=.5)
ghmma <- ricker_hmm_TMB(data=cowdat, tv.par='a',dirichlet_prior=dirpr)
ghmmb <- ricker_hmm_TMB(data=cowdat, tv.par='b',dirichlet_prior=dirpr)
ghmm <- ricker_hmm_TMB(data=cowdat, tv.par='both',dirichlet_prior=dirpr)

cp<-c(gsimple$conv_problem, 
gac$conv_problem, 
gtva$conv_problem, 
gtvb$conv_problem, 
gtvab$conv_problem, 
ghmma$conv_problem, 
ghmmb$conv_problem, 
ghmm$conv_problem) 

#-------------------------------------------------
#predictions
x_new=seq(0,max(cowdat$S)*1.5,length.out=200)
by_q=round(quantile(cowdat$by,seq(0,1,by=0.1)))
#static
pred_simple_cow=data.frame(pred=exp(gsimple$alpha-gsimple$beta*x_new)*x_new,
                       x_new=x_new,
                       gr=rep(.1,length(x_new)),
                       model="simple")
pred_ar1_cow=data.frame(pred=exp(gac$alpha-gac$beta*x_new)*x_new,
                    x_new=x_new,
                    gr=rep(.1,length(x_new)),
                    model="ar1")

rsrwa=list()
for(n in 1:length(by_q)){
  rsrwa[[n]]=exp(gtva$alpha[match(by_q[n],cowdat$by)]-gtva$beta*x_new)*x_new
}
pred_rwa_cow=data.frame(pred=do.call(c,rsrwa),
                    x_new=x_new,
                    gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                    model="rwa")
alpha_ptva=data.frame(by=cowdat$by,value=gtva$alpha,param="alpha")

rsrwb=list()
for(n in 1:length(by_q)){
  rsrwb[[n]]=exp(gtvb$alpha-gtvb$beta[match(by_q[n],cowdat$by)]*x_new)*x_new
}
beta_ptvb=data.frame(by=cowdat$by,value=gtvb$beta,param="beta")
pred_rwb_cow=data.frame(pred=do.call(c,rsrwb),
                    x_new=x_new, 
                    gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                    model="rwb")


rsrwab=list()
for(n in 1:length(by_q)){
  rsrwab[[n]]=exp(gtvab$alpha[match(by_q[n],cowdat$by)]-gtvab$beta[match(by_q[n],cowdat$by)]*x_new)*x_new
}
beta_ptvab=data.frame(by=cowdat$by,value=c(gtvab$alpha,gtvab$beta),param=rep(c("alpha","beta"),each=length(hardat$by)))

pred_rwab_cow=data.frame(pred=do.call(c,rsrwab),
                     x_new=x_new,
                     gr=rep(seq(0,1,by=0.1),each=length(x_new)),
                     model="rwab"
)


rshmma1=exp(ghmma$alpha[1]-ghmma$beta*x_new)*x_new
rshmma2=exp(ghmma$alpha[2]-ghmma$beta*x_new)*x_new
pred_hmma_cow=data.frame(pred=c(rshmma1,rshmma2),
                     x_new=x_new,
                     gr=rep(c(.1,1),each=length(x_new)),
                     model="HMMa")

gamma_df_hmma=data.frame(by=cowdat$by,gamma=ghmma$probregime[2,])

rshmmb1=exp(ghmmb$alpha-ghmmb$beta[1]*x_new)*x_new
rshmmb2=exp(ghmmb$alpha-ghmmb$beta[2]*x_new)*x_new
pred_hmmb_cow=data.frame(pred=c(rshmmb1,rshmmb2),
                     x_new=x_new,
                     gr=rep(c(1,.1),each=length(x_new)),
                     model="HMMb")

gamma_df_hmmb=data.frame(by=cowdat$by,gamma=ghmmb$probregime[2,])



rshmmab1=exp(ghmm$alpha[1]-ghmm$beta[1]*x_new)*x_new
rshmmab2=exp(ghmm$alpha[2]-ghmm$beta[2]*x_new)*x_new
pred_hmmab_cow=data.frame(pred=c(rshmmab1,rshmmab2),
                      x_new=x_new,
                      gr=rep(c(.1,1),each=length(x_new)),
                      model="HMMab")

gamma_df_hmmab=data.frame(by=cowdat$by,gamma=ghmm$probregime[2,])

pred_cow<-rbind(pred_simple_cow,
  pred_ar1_cow,
  pred_rwa_cow,
  pred_rwb_cow,
  pred_rwab_cow,
  pred_hmma_cow,
  pred_hmmb_cow,
  pred_hmmab_cow)

#--------------------------------------------
#prediction plot
pred_cow$model<-factor(pred_cow$model, levels=c("rwa", "rwb", "rwab","HMMa","HMMb", "HMMab","simple", "ar1") )

p<-ggplot(pred_cow)+
  geom_line(aes(x=x_new,y=pred,colour = as.factor(gr)),linewidth=1.5)+
  geom_point(data=cowdat, aes(x=S, y=R,alpha=by),size=1.5) +
  scale_colour_viridis_d(name='time')+
  facet_wrap(~model)+
  xlab("Spawners") + 
  ylab("Recruits")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
p
ggsave("tex/predcow.png",plot=p,width = 10, height = 7)



pred_cow_hmma<-pred_cow[pred_cow$model=="HMMa",]
library(relayer)
library(dplyr)

ggplot(pred_cow_hmma)+
  geom_line(aes(x=x_new,y=pred,colour1 = as.factor(gr)),linewidth=1.5)%>%
     rename_geom_aes(new_aes = c("colour" = "colour1")) +
  geom_point(data=cowdat, aes(x=S, y=R,colour2=by),size=2) %>%
     rename_geom_aes(new_aes = c("colour" = "colour2")) +
  scale_colour_viridis_d(aesthetics = "colour1", guide = "legend", name = "regime", option = "A",end=.8) +
  scale_colour_viridis_c(aesthetics = "colour2", guide = "legend", name = "Year", option = "D") +
  xlab("Spawners") + 
  ylab("Recruits")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(face="bold", size=14),axis.title = element_text(face="bold", size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=14),
        legend.position = "none")

x_new=seq(0,max(cowdat$S)*1.5,length.out=200)
by=cowdat$by
rshmma1=exp(ghmma$alpha[1]-ghmma$beta*x_new)*x_new
rshmma2=exp(ghmma$alpha[2]-ghmma$beta*x_new)*x_new
pred_hmma_cow_plot=data.frame(pred=c(rshmma1,rshmma2),
                     x_new=x_new,
                     by=rep(cowdat$by,each=length(x_new)),
                     model="HMMa",
                     S=rep(cowdat$S,each=length(x_new)),
                     R=rep(cowdat$R,each=length(x_new)))


pred_cow_hmma<-pred_cow[pred_cow$model=="HMMa",]
head(pred_cow_hmma)
ggplot(pred_hmma_cow_plot)+
  geom_line(aes(x=x_new,y=pred,alpha = as.factor(by)),linewidth=1.5)+
  geom_point(aes(x=S, y=R,colour=by),size=1.5) +
  #scale_colour_viridis_d(name='regime')+
  scale_colour_viridis_d(name='Year')+
  xlab("Spawners") + 
  ylab("Recruits")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))



#==========================================
#parameter trends

dfpar_cow<-data.frame(by=rep(cowdat$by,16),
                  value=c(
                    rep(gsimple$alpha,length(cowdat$by)),
                    rep(gac$alpha,length(cowdat$by)),
                    gtva$alpha,
                    rep(gtvb$alpha,length(cowdat$by)),
                    gtvab$alpha,
                    ghmma$alpha[ghmma$regime],
                    rep(ghmmb$alpha,length(cowdat$by)),
                    ghmm$alpha[ghmm$regime],
                    1/rep(gsimple$beta,length(cowdat$by)),
                    1/rep(gac$beta,length(cowdat$by)),
                          1/rep(gtva$beta,length(cowdat$by)),
                          1/gtvb$beta,
                          1/gtvab$beta,
                          1/rep(ghmma$beta,length(cowdat$by)),
                          1/ghmmb$beta[ghmmb$regime],
                          1/ghmm$beta[ghmm$regime]),
                    param=rep(c("log_a","Smax"), each=length(cowdat$by)*8),
                    model=factor(rep(rep(c("simple","AR1","RWa","RWb","RWab",
                                "Shifta","Shiftb","Shiftab"),2),each=length(cowdat$by)),
                                levels=c("simple","AR1","RWa","RWb","RWab",
                                         "Shifta","Shiftb","Shiftab")))




pp<-ggplot(dfpar_cow) +   
  geom_line(aes(x=by,y=value,color=by ),
            linewidth=1.4) +  
  mytheme+ 
  ylab("parameters") +
  xlab("year") +
  facet_grid(param~model, scales="free_y")+
  scale_color_viridis_c(begin=.1, end=.8, guide="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
pp

ggsave("tex/paramcow.png",plot=pp,width = 12, height = 5)


dfpar_cow_hmma<-dfpar_cow[dfpar_cow$model=="Shifta"&dfpar_cow$param=="log_a",]

ggplot(dfpar_cow_hmma) +   
  geom_line(aes(x=by,y=value ),
            linewidth=1.4) +  
  geom_point(aes(x=by,y=value,color=by ),size=3            ) + 
  mytheme+ 
  ylab("log(Alpha)") +
  xlab("Year") +
  scale_color_viridis_c()+
  theme(legend.position="right")
#model selection criteria
#lfo
#lfo
    glfostatic <- tmb_mod_lfo_cv(data=cowdat,model='static', L=15)
    glfoac <- tmb_mod_lfo_cv(data=cowdat,model='staticAC', L=15)
    glfoalpha <- tmb_mod_lfo_cv(data=cowdat,model='rw_a', siglfo="obs", L=15)
    glfobeta <- tmb_mod_lfo_cv(data=cowdat,model='rw_b', siglfo="obs", L=15)
    glfoalphabeta <- tmb_mod_lfo_cv(data=cowdat,model='rw_both', siglfo="obs", L=15)
    glfohmma <- tmb_mod_lfo_cv(data=cowdat,model='HMM_a', L=15)
    glfohmmb <- tmb_mod_lfo_cv(data=cowdat,model='HMM_b', L=15)
    glfohmm <- tmb_mod_lfo_cv(data=cowdat,model='HMM', L=15)
 lfores<-c(sum(glfostatic$lastparam), 
        sum(glfoac$lastparam), 
        sum(glfoalpha$lastparam), 
        sum(glfobeta$lastparam), 
        sum(glfoalphabeta$lastparam), 
        sum(glfohmma$lastregime_pick),
        sum(glfohmmb$lastregime_pick),
        sum(glfohmm$lastregime_pick))
                           


df<-data.frame(model=c("simple","AR1","RW a","RW b","RW a and b",
                    "Shift a","Shift b","Shift a and b"),
           AIC=round(c(gsimple$AICc, gac$AICc,
                 gtva$AICc,gtvb$AICc, gtvab$AICc,
                 ghmma$AICc,ghmmb$AICc, ghmm$AICc)-min(
                c(gsimple$AICc, gac$AICc,
                     gtva$AICc,gtvb$AICc, gtvab$AICc,
                 ghmma$AICc,ghmmb$AI, gtvab$AICc)),2),
           BIC=round(c(gsimple$BIC, gac$BIC,
                 gtva$BIC,gtvb$BIC, gtvab$BIC,
                 ghmma$BIC, ghmmb$BIC, ghmm$BIC)-
             min(c(gsimple$BIC, gac$BIC,
                   gtva$BIC,gtvb$BIC, gtvab$BIC,
                   ghmma$BIC, ghmmb$BIC, ghmm$BIC)),2),
             LFO=round(c(sum(glfostatic$lastparam), 
                       sum(glfoac$lastparam), 
                       sum(glfoalpha$lastparam), 
                       sum(glfobeta$lastparam), 
                       sum(glfoalphabeta$lastparam), 
                       sum(glfohmma$lastregime_pick),
                       sum(glfohmmb$lastregime_pick),
                       sum(glfohmm$lastregime_pick)),2))
             

df$AIC = ifelse(df$AIC ==0, paste0("\\textbf{", df$AIC, "}"), df$AIC)
df$BIC = ifelse(df$BIC ==0, paste0("\\textbf{", df$BIC, "}"), df$BIC)


df$AIC[cp]<-NA
df$BIC[cp]<-NA

df$LFO = ifelse(df$LFO ==max(df$LFO, na.rm=T), paste0("\\textbf{", df$LFO, "}"), df$LFO)



print(xtable(df,caption="Parameter trajectories for stationary 
  and time varying models considered for Cowichan Chinook. Minimum 
  AIC and BIC value indicate preferred model, while maximum LFO indicate 
  preferred model. Missing values indicate lack of model convergence",
label="cowcrit"), caption.placement="top", file="tex/cowcrit.tex",
sanitize.text.function = identity,NA.string = "-")




