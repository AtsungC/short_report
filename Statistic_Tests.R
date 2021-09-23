library(tidyverse)


################## trial #################
#### variance and t tests with df_trial (6,8 wk old)####
# raw data = actual figure result of the tests
# sign data = sign indicate the result of the tests (alpha = 0.05) 
# ft_test(data_frame,TX_group1,TX_group2)

unique(df_trial$Group) #"HOM_Untreated"   "WT_Tanganil_8wk"  "WT_Untreated"     "HOM_Tanganil_8wk"
trial_test_68wk <- ft_test(df_trial,"HOM_Untreated","WT_Untreated")
trial_test_68wk <- ft_test(df_trial,"WT_Tanganil_8wk","HOM_Tanganil_8wk")

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_test_68wk, "trial_t_test_68wk_raw.csv", row.names=FALSE, quote=FALSE) 

# remove the variance test result
# the the table of (sign/not sign)
trial_test_68wk <- trial_test_68wk %>% select(-2,-4)
trial_test_68wk[,2:3] <- sapply(trial_test_68wk[,2:3],as.numeric) 
trial_test_68wk[,-1][trial_test_68wk[,-1]<=0.05] <- 'sign'
trial_test_68wk[,-1][trial_test_68wk[,-1]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_test_68wk, "trial_t_test_68wk_sig.csv", row.names=FALSE, quote=FALSE) 


#### aov, repeated aov and trend : df_trial####
# data = df_trial
# raw data = actual figure result of the tests
# sign data = sign indicate the result of the tests (alpha = 0.05) 
# three weeks 6,7,8
#aov_test (wt_group,npc_group,data)

trial_repeated_aov <- aov_test("WT_Untreated","HOM_Untreated",df_trial)
trial_repeated_aov <- aov_test("WT_Tanganil_8wk","HOM_Tanganil_8wk",df_trial)

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_repeated_aov, "trial_aov_test_raw.csv", row.names=FALSE, quote=FALSE) 
trial_repeated_aov[,2:5][trial_repeated_aov[,2:5]<=0.05] <- 'sign'
trial_repeated_aov[,2:5][trial_repeated_aov[,2:5]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_repeated_aov, "trial_aov_test_sign.csv", row.names=FALSE, quote=FALSE)

# combine the result with variance&t tests (sign version)
# summary table
identical(trial_test_68wk[[1]],trial_repeated_aov[[1]]) # true
test_result <- full_join(trial_test_68wk,trial_repeated_aov,by=c('variables'='variable'))
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(test_result, "trial_test_result.csv", row.names=FALSE, quote=FALSE)


## trend analysis ##
# linear and quadratic regression
#trend_text (wt, npc, data frame)
var_selected <- trial_repeated_aov[[1]]


# 

trial_mixed_aov <- aov_mixed("WT_Untreated","HOM_Untreated",df_trial) %>% na.omit()
# trial_mixed_aov[,2:4][trial_mixed_aov[,2:4]<=0.10] <- 'sign' # increase type I error to 10% to increase the significance
# trial_mixed_aov[,2:4][trial_mixed_aov[,2:4]!='sign'] <- 'not sign'
trial_mixed_aov <- gather(trial_mixed_aov,ana,p_value,inter_p:Group_p,factor_key = T) 
var_selected <- trial_mixed_aov %>% filter(ana=='inter_p',p_value<=0.1)
trial_mixed_aov <- trial_mixed_aov %>% filter(variable %in% var_selected[[1]],p_value<=0.1) %>% arrange(p_value)
trial_mixed_aov %>% ggplot(aes(x=ana,y=reorder(variable,p_value),fill=p_value))+
  geom_tile()+
  scale_fill_viridis_c(option = "B", direction = -1)+
  labs(title='Mixed ANOVA model',subtitle = 'WT vs NPC (From 6 to 8wk)',fill='P value')+
  theme(axis.title.y =element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 5))

trial_trend_aov <- trend_test("WT_Untreated","HOM_Untreated",df_trial)
# trial_trend_aov <- trend_test("WT_Tanganil_8wk","HOM_Tanganil_8wk",df_trial)  
trial_trend_aov<- gather(trial_trend_aov,ana,p_value,WT_linear_p:NPC_quad_p,factor_key = T) 
trial_trend_aov <- trial_trend_aov %>% filter(abs(p_value)<=0.1)
trial_p <- rbind(trial_trend_aov,trial_mixed_aov)

P <- trial_p %>%  mutate(name = fct_relevel(ana,'inter_p','Age_p','Group_p','WT_linear_p','NPC_linear_p','WT_quad_p','NPC_quad_p')) %>%
  ggplot(aes(x=name,y=reorder(variable,p_value),fill=p_value))+
  geom_tile()+
  theme_classic()  +
  scale_fill_viridis_c( direction = -1)+
  labs(title='Mixed ANOVA model',subtitle = 'WT vs NPC (From 6 to 8wk)',fill='P value')+
  theme(axis.title.y =element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 8)) +
  NULL

  ggsave(P,filename='heatmap_p.png')

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_trend_aov, "trial_trend_raw.csv", row.names=FALSE, quote=FALSE) 

#### linear regression to prediciton at 9wk with CI ####
trial_p %>% filter(ana == 'NPC_linear_p')
summary(trial_p)
linear_var <- unique(as.character(trial_p$variable[trial_p$ana=='NPC_linear_p']))

p9 <- pred_9wk(df_trial,'HOM_Untreated',linear_var) %>% mutate(Group='pred') %>% select("Group","variable","pred", "se","ci_l","ci_h",-'Age')

fixed_9 <- df_trial  %>% select(ID,Age,Group,linear_var)%>% filter(Group%in% c("WT_Untreated","HOM_Tanganil_8wk"),Age==9)
fixed_9<- gather(fixed_9,variable,fixed,linear_var,factor_key = T) %>% select(variable,fixed,Group)

fixed_9 <- fixed_9 %>% group_by(Group,variable) %>% summarise(pred=mean(fixed),se=sd(fixed),ci_l=pred-1.96*se,ci_h=pred+1.96*se)
length(unique(fixed_9$variable))
pd_9 <- rbind(fixed_9,p9)
# pd_9 <- pd_9 %>% mutate(type_v=ifelse(str_detect(pd_9$variable,'%'),1,0))
pd_9$type_v <-  as.factor(c('percent','unit')[as.numeric(str_detect(pd_9$variable, '%'))+1])

pd_9_ctr <- pd_9 %>% filter(Group %in% c('HOM_Tanganil_8wk','WT_Untreated')) %>% select('Group','variable','pred') %>% mutate(ci_l=pred,ci_h=pred)
pd_9_tx <- pd_9 %>% filter(Group =='pred') %>% select('Group','variable','pred','ci_l','ci_h') 
pd_9 <- rbind(pd_9_ctr,pd_9_tx)
pd_9 <- pd_9 %>% mutate(tx_ben=ifelse(pd_9$ci_l[pd_9$Group=='pred']< pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk'] & pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk']<pd_9$ci_h[pd_9$Group=='pred'],1,0))

pd_9 <- pd_9 %>% mutate(tx_ben=ifelse(pd_9$pred[pd_9$Group=='WT_Untreated']< pd_9$pred[pd_9$Group=='pred'],ifelse(pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk']<pd_9$pred[pd_9$Group=='pred']&pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk']>pd_9$pred[pd_9$Group=='WT_Untreated'],1,0),ifelse(pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk']>pd_9$pred[pd_9$Group=='pred']&pd_9$pred[pd_9$Group=='HOM_Tanganil_8wk']<pd_9$pred[pd_9$Group=='WT_Untreated'],1,0)))
#### check the order with or without the log transformation and the magnitude ####
pd_9 <- pd_9[!str_detect(pd_9$variable,'Intensity'),]
pd_9_h <- pd_9[pd_9$Group=='pred'&pd_9$pred>=10,]$variable
pd_9_l <- pd_9[pd_9$Group=='pred'&pd_9$pred<=10,]$variable
# pd_9 <- pd_9 %>% filter(variable!='StepSequence_RB_(%)')
library(skimr)
pd_9 %>% group_by(Group,variable) %>% skim_to_wide()

pd_9  %>% filter(tx_ben==1 & variable %in% pd_9_h) %>%ggplot(aes(y=variable,x=pred,color=Group))+
  geom_point()+
  geom_errorbar(aes(xmin=ci_l,xmax=ci_h),width=0)+ 
  NULL
pd_9  %>% filter(tx_ben==1 & variable %in% pd_9_l) %>%ggplot(aes(y=variable,x=pred,color=Group))+
  geom_point()+
  geom_errorbar(aes(xmin=ci_l,xmax=ci_h),width=0) + 
  NULL
pd_9 %>% filter(tx_ben==1)  %>%ggplot(aes(y=variable,x=pred,color=Group))+
  geom_point()+
  geom_errorbar(aes(xmin=ci_l,xmax=ci_h),width=0) + 
  NULL
p1 <- pd_9 %>% filter(type_v=='percent') %>%
  ggplot(aes(y=variable,x=pred,color=Group))+
  geom_point()+
  geom_errorbar(aes(xmin=ci_l,xmax=ci_h), position='dodge2') + 
  NULL
p2 <- pd_9 %>% filter(type_v=='unit') %>%
  ggplot(aes(y=variable,x=pred,color=Group))+
  geom_point(position = position_dodge2(width=0.9), 
             aes(y=variable, color=Group)) +
  geom_errorbar(aes(y=variable,xmin=ci_l,xmax=ci_h),
                position=position_dodge2(width=0.9)) + 
  NULL
  # coord_cartesian(ylim=c(-12,12))
  print(p)
  #ggsave(p,filename = 'pred vs fixed.png')
#
  pred_list %>% ggplot(aes(y=variable,x=pred))+
    geom_point(position = position_dodge2(width=0.9), 
               aes(y=variable)) +
    geom_errorbar(aes(y=variable,xmin=ci_l,xmax=ci_h),
                  position=position_dodge2(width=0.9)) + 
    NULL
  
  

#
trial_trend_aov[,-1][trial_trend_aov[,-1]<=0.05] <- 'sign'
trial_trend_aov[,-1][trial_trend_aov[,-1]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(trial_trend_aov, "trial_trend_sign.csv", row.names=FALSE, quote=FALSE) 

# combine the result with previous summary table 
identical(test_result[[1]],trial_trend_aov[[1]]) # true
test_result_1 <- full_join(test_result,trial_trend_aov,by=c('variables'='variable'))
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(test_result_1, "trial_test_result.csv", row.names=FALSE, quote=FALSE)




################## run ###################
#### variance and t tests with df_run (6,8 wk old)####

run_test_group <- list()
var_selected <- colnames(df_run)
run_test_group <- ft_test(df_trial,"HOM_Untreated","WT_Untreated")
run_test_group <- ft_test(df_trial,"WT_Tanganil_8wk","HOM_Tanganil_8wk")

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_test_68wk, "t_test_68wk_raw.csv", row.names=FALSE, quote=FALSE) 
# remove the variance test result
run_test_68wk <- run_test_68wk %>% select(-2,-4)
run_test_68wk[,2:3] <- sapply(run_test_68wk[,2:3],as.numeric)
run_test_68wk[,-1][run_test_68wk[,-1]<=0.05] <- 'sign'
run_test_68wk[,-1][run_test_68wk[,-1]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_test_68wk, "t_test_68wk_sig.csv", row.names=FALSE, quote=FALSE) 

#### matched t-test ####
da1 <- da %>% select(ID,Age,Group,Trial,Run,`Run_Duration_(s)`) %>% spread(Age,`Run_Duration_(s)`)
t.test(as.vector(na.omit(da1$`6`[da1$ID=='386.4a'])),as.vector(na.omit(da1$`8`[da1$ID=='386.4a'])),paired = T,var.equal = F,alternative = 'two.sided',mu=0)

####aov, repeated aov and trend : df_run####
#df_run # all data in run_statistic
var_selected <- run_test_68wk[[1]]# consistance with t-test

run_repeated_aov <- aov_test("WT_Untreated","HOM_Untreated",df_run)
run_repeated_aov <- aov_test("WT_Tanganil_8wk","HOM_Tanganil_8wk",df_run)

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_repeated_aov, "run_aov_test_raw_tg.csv", row.names=FALSE, quote=FALSE) 

run_repeated_aov[,2:5][run_repeated_aov[,2:5]<=0.05] <- 'sign'
run_repeated_aov[,2:5][run_repeated_aov[,2:5]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_repeated_aov, "run_aov_test_sign_tg.csv", row.names=FALSE, quote=FALSE) 
identical(run_test_68wk[[1]],run_repeated_aov[[1]]) # true

test_result <- full_join(run_test_68wk,run_repeated_aov,by=c('variables'='variable'))
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(test_result, "run_test_result_tg.csv", row.names=FALSE, quote=FALSE)

## trend analysis ##
# linear and quadratic regression
#trend_text (wt, npc, data frame)
var_selected <- run_repeated_aov[[1]]
run_trend_aov <- trend_test("WT_Untreated","HOM_Untreated",df_run)
run_trend_aov <- trend_test("WT_Tanganil_8wk","HOM_Tanganil_8wk",df_run)

setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_trend_aov, "run_trend_raw_tg.csv", row.names=FALSE, quote=FALSE) 

run_trend_aov[,-1][run_trend_aov[,-1]<=0.05] <- 'sign'
run_trend_aov[,-1][run_trend_aov[,-1]!='sign'] <- 'not sign'
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(run_trend_aov, "run_trend_sign_tg.csv", row.names=FALSE, quote=FALSE) 

identical(test_result[[1]],run_trend_aov[[1]]) # true
test_result_1 <- full_join(test_result,run_trend_aov,by=c('variables'='variable'))
setwd('/Volumes/KINGSTON/Ecem_Catwalk/')
write.csv(test_result_1, "run_test_result_tg.csv", row.names=FALSE, quote=FALSE)



