library(tidyverse)



library(lme4)
library(nlme)

aov_mixed <- function(G1_WT,G2_NPC,df){

# df <- df_trial
# G1_WT <- "WT_Untreated"
# G2_NPC <- "HOM_Untreated"
  trial_repeated_aov <- list()
  df_aov <- df %>% filter(Group%in% c(G1_WT,G2_NPC),Age %in% c(6,7,8))
  var_selected <- colnames(df_aov)
  for (i in 13:length(var_selected)) {
    tryCatch({
      #i<- 258    #gives StepSequence_RegularityIndex
      j <- i-12
      df <- df_aov%>% select('ID','Age','Group',var_selected[i])
      names(df) <- c('ID','Age','Group','Y')
      if(all(is.na(df[[4]]))| length(unique(df[[4]]))==1){
        trial_repeated_aov[[j]] <- data.frame(variable =var_selected[i],
                                              inter_p  =NA,
                                              Age_p    =NA,
                                              Group_p  =NA
                                              )
        next}
      # df$Age   <- factor(df$Age)
      df$Group <- factor(df$Group)
      df$ID    <- as.factor(df$ID)
      df$AgeOrd<- as.ordered(df$Age) ### oordered factors allow you to fit linear and quadratic
      ### effects - this might be important later but will decrease the significance
      ### with respect to that from the model including Age as a numeric variable
      df<- nlme::groupedData(form= Y ~ 1 | ID, data=df)
      #normality <- shapiro.test(df$Y)$p.value
      ### variances across groups should be homogeneous - check with residuals
      aov_re <- aov(data=df,Y~Age*Group+Error(ID)) # two-way + repeated anova = mixed model
      # lm.fit<- lm( Y ~ Age*Group, data =df) ## model without random effects on ID
      # lme.fit<- nlme::lme(fixed=Y ~ Age*Group, random=~1, data=df, method='ML') # fit with ML, not REML, to compare BIC
      # lme.resid<-  residuals(lme.fit, level=0:1) 
      # ### gives an n x 2 matrix with cols fixed, ID,
      # ### the fixed residuals are the population (fixed effects residuals), the ID one
      # ### is the indiviual residuals for the random effects within ID
      # n.obs<- length(lme.resid)/2 ### number of rows in residual object
      # df.resid<- data.frame(residual=as.vector(lme.resid), 
      #                          residual.type=as.factor(rep(1:2, rep(n.obs,2))))
      # attr(df.resid$residual.type, 'levels')<- c('fixed','ID')
      # ### plot residuals 
      # df.resid %>% ggplot(aes(sample=residual)) +
      #              stat_qq() +
      #              stat_qq_line(col='red') + 
      #              facet_wrap(~residual.type, ncol=3) + 
      #              ggtitle(var_selected[i]) + 
      #              NULL
      # 
      inter_p <- unlist(summary(aov_re)[[2]])[[14]]
      Age_p   <- unlist(summary(aov_re)[[2]])[[13]]
      Group_p <- unlist(summary(aov_re)[[1]])[[9]]
      res_aov <- residuals(aov_re)
      # lme_re <- nlme::lme(colnames(df)[[3]]~colnames(df)[[2]], random = ~ 1 | colnames(df)[[1]], data=df)
      # lme_re <- nlme::lme(`RF_StandIndex_Mean` ~ Age, random = ~ 1 | ID/Age, data=df)
      
      # library(multcomp)
      # sum_lme <- summary(multcomp::glht(lme_re, linfct=mcp(Age="Tukey")), test = adjusted(type = "bonferroni"))

      trial_repeated_aov[[j]] <- data.frame(variable=var_selected[i],inter_p=round(inter_p ,3),Age_p=round(Age_p ,3),Group_p=round(Group_p ,3))
    },error=function(e){cat("ERROR :",conditionMessage(e),'variable_i: ',var_selected[[i]],"_",i, "\n")})
  } 
  trial_repeated_aov <- do.call(rbind.data.frame,trial_repeated_aov) 
  trial_repeated_aov[,2:3] <- sapply(trial_repeated_aov[,2:3],as.numeric) 
  trial_repeated_aov <- arrange(trial_repeated_aov,inter_p,Age_p,Group_p)
  return(trial_repeated_aov)
}

# saveRDS(object ='aov_mixed',file = 'aov_mixed.rds')

trend_test <- function(G1_WT,G2_NPC,df){
  library(nlme)
  # c1 <- c(-1,0,1)
  # c2 <- c(0.5,-1,0.5)
  # mat <- cbind(c1,c2)
  var_selected <- trial_mixed_aov[[1]]
  # df <- df_trial
  # G1_WT <- "WT_Untreated"
  # G2_NPC <- "HOM_Untreated"
  summary(df$Age)
  # 
  df_aov <- df %>% filter(Group%in% c(G1_WT,G2_NPC),Age %in% c(6,7,8))
  summary(df_aov$Age)
  trial_trend_aov <- list()
  for (i in 1:length(var_selected)) {
    tryCatch({
      j <- i
      df <- df_aov%>% filter(Group==G1_WT) %>% select('ID','Age',var_selected[i])
      if(all(is.na(df[[3]]))){
        trial_trend_aov[[j]] <- data.frame(variable=var_selected[i],WT_linear_p=NA,WT_quad_p=NA,
                                           NPC_linear_p=NA,NPC_quad_p=NA)
        next}
      df$Age <- as.ordered(df$Age) # generate the contrast 
      summary(df$Age)
      class(df$Age)
      # contrasts(df$Age) <- mat
      names(df)[3]<- "Y"
      df$ID<- as.factor(df$ID)
      aov_rep <- aov(data=df,Y~Age+Error(ID))
   
      df<- nlme::groupedData(Y ~ Age | ID, data=df) #declares a nested structure within ID
      lme_rep<- nlme::lme(Y~Age, random=~1, data=df) #note that the p-values are equal to those from aov_rep
      L_s<- sign(lme_rep$coefficients$fixed[[2]])
      Q_s<- sign(lme_rep$coefficients$fixed[[3]])
      
    
      p_value_linear <- summary(aov_rep,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[2]][[1]][[2,5]]*L_s
      p_value_quad <- summary(aov_rep,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[2]][[1]][[3,5]]*Q_s
      # aov <- aov(data=df,df[[3]]~Age)
      # p_value_linear <- summary(aov,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[1]][[2,5]]
      # p_value_quad <- summary(aov,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[1]][[3,5]]
      
     
     
      # df<- nlme::groupedData(Y ~ Age | ID, data=df) #declares a nested structure within ID
      # lme_rep<- nlme::lme(Y~Age, random=~1, data=df) #note that the p-values are equal to those from aov_rep
      # beta_L_s<- sign(lme_rep$coefficients$fixed[2])
      # beta_Q_s<- sign(lme_rep$coefficients$fixed[3])

      
      df_npc <- df_aov%>% filter(Group==G2_NPC) %>% select('ID','Age',var_selected[i])
      df_npc$Age <- as.ordered(df_npc$Age)
      # contrasts(df_npc$Age) <- mat
      names(df_npc)[3] <- 'Y'
      aov_rep <- aov(data=df_npc,Y~Age+Error(ID))
     
      df_npc<- nlme::groupedData(Y ~ Age | ID, data=df) #declares a nested structure within ID
      lme_rep<- nlme::lme(Y~Age, random=~1, data=df) #note that the p-values are equal to those from aov_rep
      L_s<- sign(lme_rep$coefficients$fixed[[2]])
      Q_s<- sign(lme_rep$coefficients$fixed[[3]])
    
      
      
      p_value_linear_n <- summary(aov_rep,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[2]][[1]][[2,5]]*L_s
      p_value_quad_n <- summary(aov_rep,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[2]][[1]][[3,5]]*Q_s
      # aov <- aov(data=df_npc,df_npc[[3]]~Age)
      # p_value_linear_n <- summary(aov,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[1]][[2,5]]
      # p_value_quad_n <- summary(aov,split=list(Age=list('Linear'=1,'Quadratic'=2)))[[1]][[3,5]]

      
      trial_trend_aov[[j]] <- data.frame(variable=var_selected[i],WT_linear_p=round(p_value_linear,3),WT_quad_p=round(p_value_quad,3),
                                         NPC_linear_p=round(p_value_linear_n,3),NPC_quad_p=round(p_value_quad_n,3))
    },error=function(e){cat("ERROR :",conditionMessage(e),'variable_i: ',var_selected[i],"_",i, "\n")})
  }
  trial_trend_aov <- do.call(rbind.data.frame,trial_trend_aov)
  trial_trend_aov[,2:5] <- sapply(trial_trend_aov[,2:5],as.numeric) 
  
  return(trial_trend_aov)
}
 
t_9 <- function(G1_WT,G2_NPC,df,variable_selected){
  
  # df <- df_trial
  # G1_WT <- "WT_Untreated"
  # G2_NPC <- "HOM_Tanganil_8wk"
  trial_repeated_aov <- list()
  df_aov <- df %>% filter(Group%in% c(G1_WT,G2_NPC),Age ==9)
  var_selected <- as.character(variable_selected)
  for (i in 1:length(var_selected)) {
    tryCatch({
      #i<- 258    #gives StepSequence_RegularityIndex
      
      df <- df_aov%>% select('ID','Group',var_selected[i])
      names(df) <- c('ID','Group','Y')
      if(all(is.na(df[[3]]))| length(unique(df[[3]]))==1){
        trial_repeated_aov[[i]] <- data.frame(variable =var_selected[i],
                                              p_value =NA,
                                              conf_l=NA,
                                              conf_h=NA,estimate_dff=NA)
        next}
      
      f <- var.test(df[,3][df[,2]==G1_WT],df[,3][df[,2]==G2_NPC],ratio = 1,alternative = 'two.sided',na.action=na.omit())
      if(!is.na(f[[3]])){
        ifelse(f[[3]]>=0.05,t <- t.test(df[,3][df[,2]==G1_WT],df[,3][df[,2]==G2_NPC],var.equal = T,alternative = 'two.sided',mu=0,na.action=na.omit()),
               t <-t.test(df[,3][df[,2]==G1_WT],df[,3][df[,2]==G2_NPC],var.equal = F,alternative = 'two.sided',mu=0,na.action=na.omit()))
      }
      
      trial_repeated_aov[[i]] <- data.frame(variable=var_selected[i],p_value=round(t$p.value ,3),conf_l=t$conf.int[1],conf_h=t$conf.int[2],estimate_dff=((t$estimate[[1]]-t$estimate[[2]])))
    },error=function(e){cat("ERROR :",conditionMessage(e),'variable_i: ',var_selected[[i]],"_",i, "\n")})
  } 
  trial_repeated_aov <- do.call(rbind.data.frame,trial_repeated_aov) 
  return(trial_repeated_aov)
}


pred_9wk <- function(df,G1,var){
  df1 <- df %>% select(ID,Age,Group,var) %>% filter(Age %in% c(6,7,8),Group==G1)
  pred_list <- list()
  for (i in 1:length(var)) {
    df2 <- df1 %>% select(ID,Age,var[[i]])
    # df2$Age <- as.ordered(df2$Age) # generate the contrast 
 
    # contrasts(df$Age) <- mat
    names(df2)[3]<- "Y"
    df2$ID<- as.factor(df2$ID)
    if(min(df2$Y)>0) df2$Y<- log(df2$Y)
    df2<- nlme::groupedData(Y ~ Age | ID, data=df2) #declares a nested structure within ID
    lme_rep<- nlme::lme(Y~Age, random=~1, data=df2) #note that the p-values are equal to those from aov_rep
    age_9 <- data.frame(Age=9)
    pd <- bernr::bolker_ci(lme_rep,age_9)
    pd <- cbind(data.frame(variable=var[[i]],pd))
    pd.new<- pd
    pd.new$pred<- exp(pd$pred); pd.new$ci_l<- exp(pd$pred-1.96*pd$se); pd.new$ci_h<-exp(pd$pred+1.96*pd$se)
    pred_list[[i]] <- pd.new
  }
  pred_list <- do.call(rbind.data.frame,pred_list)
}




