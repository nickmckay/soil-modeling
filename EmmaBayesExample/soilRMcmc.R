# Bayesian radiocarbon modeling -------------------------------------------
library(dplyr)
library(SoilR)
library(BayesianTools)
library(foreach)
library(ggplot2)
library(tidyr)
rm(list=ls()) # clear the workspace

#Set plot themes 
ts_theme <- theme_classic() +
  theme(plot.title = element_text(size = 26, hjust = 0.5),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.key = element_rect(fill = "transparent"),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks=element_line(color = "black", linewidth=1),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.background = element_rect(color="black", fill="grey83", size=1.5, linetype="solid"))
theme_set(ts_theme)

treatPal <- c("#808080","#377EB8", "#E41A1C") 
depthPal <- c("#9AA582FF", "#657359FF", "#8B775FFF", "#6B472CFF")

# File paths --------------------------------------------------------------
MCMC_path <- '/Users/emmalathrop/Documents/NAU/SchuurLab/RCode/SoilR/MCMCOutput/'
manDir <- "/Users/emmalathrop/Documents/NAU/SchuurLab/Manuscripts/Chapter1/figures/RFigs/"

# Import soil measurements ---------------------------------------------------------

c14Dat <- read.csv('/Volumes/Schuur Lab/2020 New_Shared_Files/DATA/CiPEHR & DryPEHR/Soil Cores/2022/Datasets/c14Summary20250507.csv') %>% 
  dplyr::select(-X) %>% 
  mutate(depth = factor(depth, levels = c("Surface organic", "Deep organic", "Mineral", "Deep mineral"),
                        labels = c("so", "do", "m", "dm")),
         treatment = ifelse(treatment == "2009", "Initial", treatment),
         treatment = ifelse(treatment == "2022: Control", "Control", treatment),
         treatment = ifelse(treatment == "2022: Warming", "Warming", treatment),
         year = ifelse(treatment == "Initial", 2009, 2022)) %>% 
  filter(depth != "dm")

cStockDat <- read.csv('/Volumes/Schuur Lab/2020 New_Shared_Files/DATA/CiPEHR & DryPEHR/Soil Cores/2022/Datasets/stocks/stocksWithMin4Depths20250507.csv') %>% 
  dplyr::select(year, treatment:stock.m) %>% 
  tidyr::pivot_longer(cols  = -c(year, treatment), names_to = c(".value", "depth"), 
                      names_sep = "\\.") %>%
  mutate(treatment = ifelse(year == 2009, "Initial", treatment),
         treatment = ifelse(treatment == "c", "Control", treatment),
         treatment = ifelse(treatment == "w", "Warming", treatment)) %>% 
  filter(depth != "dm") %>% 
  group_by(year, treatment, depth) %>%
  summarise(meanCStock = mean(stock, na.rm = TRUE),
            sdCStock = sd(stock, na.rm = TRUE),
            seCStock = sdCStock/sqrt(length(!is.na(stock)))) %>% 
  ungroup()

datComb <- left_join(cStockDat, c14Dat, by = c("year", "treatment", "depth")) %>% 
  dplyr::select(-seCStock, -se14c) %>% 
  dplyr::rename(layer = depth, C14 = mean14c, C_14_sd = sd14c, C_stock = meanCStock, C_stock_sd = sdCStock) %>%
  mutate(layer = ifelse(layer == "m", "min", layer),
         layer = factor(layer, levels = c("so", "do", "min"))) %>% 
  arrange(year, treatment, layer) %>% 
  dplyr::mutate(layer = as.character(layer))

# Create a vector of atmospheric 14C values -------------------------------

#Forecast from 2009 to 2022 to get atmospheric 14C values
atmCip <- read.csv('/Volumes/Schuur Lab/2020 New_Shared_Files/DATA/Gradient/Radiocarbon/515_EML_Radiocarbon_Atm_Up to Date.csv', skip = 1) %>% 
  dplyr::select(yearly.average:X.2) %>% 
  dplyr::rename(Atmosphere = X.2, YEAR = yearly.average) %>% 
  filter(YEAR > 2008 & YEAR < 2023)

atmIn <- rbind(subset(C14Atm_NH, YEAR < 2009),atmCip)
years=seq(2009,2022)

# Model pool change since 2009 ------------------------------------
#source code: https://rdrr.io/rforge/SoilR/src/R/ThreepSeriesModel14.R
#general model: https://rdrr.io/rforge/SoilR/man/94.html#heading-0
get_likelihood=function(pars){
  
  model_output=ThreepSeriesModel14(
    t=years,ks=c(k1=pars[1], k2=pars[2], k3=pars[3]),
    C0=dat$C_stock[dat$year==2009], F0_Delta14C=dat$C14[dat$year==2009],
    In=pars[4], a21=pars[5], a32=pars[6],
    inputFc=atmIn[atmIn$YEAR%in%2009:2022,],pass=T
  )
  
  C14=getF14(model_output)[years%in%c(2022),]
  C_stock=getC(model_output)[years%in%c(2022),]
  
  ### change this to merge properly with model output
  C14_residuals=C14-dat$C14[dat$year%in%c(2022)]
  C_stock_residuals=C_stock-dat$C_stock[dat$year%in%c(2022)]
  
  ### replace sd with observed sd of 14C and C_stock
  C14_likelihood=sum(dnorm(C14_residuals,0,dat$C_14_sd*0.4,log=T))
  C_stock_likelihood=sum(dnorm(C_stock_residuals,0,dat$C_stock_sd*0.4,log=T))
  
  total_likelihood=sum(C14_likelihood,
                       C_stock_likelihood,
                       na.rm=T)
  
  return(total_likelihood)
}



# Run MCMC Control ----------------------------------------------------------------
prior=createUniformPrior(lower=c(0.002,0.0006,0.0002,
                                 10, 0.0025,0.0001),
                         upper=c(0.2,0.06,0.06,
                                 205,0.5,0.1))
settings = list(iterations = 5000 ,message=T)
bayesianSetup <- createBayesianSetup(get_likelihood, 
                                     prior=prior, 
                                     names = c('k1','k2','k3','input', 'a21','a23')
)
dat <-  data.frame(subset(datComb, year == 2009 | treatment == "Control") %>% 
                     dplyr::select(-treatment) %>% 
                     mutate(C_stock = C_stock*1000,
                            C_stock_sd = C_stock_sd*1000))
#test this function works
get_likelihood(c(0.1,0.01,0.01,200,0.1,0.01))

#start running model
out <- runMCMC(bayesianSetup = bayesianSetup,
               sampler = "DREAMzs",
               settings = settings)
saveRDS(out, file=paste0(MCMC_path,'mcmc_Out_Control',format(Sys.Date(), "%Y%m%d%S"), Sys.time(),'.RDS'))

#out <- readRDS("~/Documents/NAU/SchuurLab/RCode/SoilR/MCMCOutput/mcmc_Out_Control20250508002025-05-08 08:16:33.759319.RDS")

plot(out)
burn_in=1500
out_sample=reshape2::melt(variable.name='variable_name',
                          data.frame(getSample(out,start=burn_in)))
MAP_params=data.frame(MAP=MAP(out,start=burn_in)$parametersMAP,
                      variable_name=rownames(data.frame(MAP(out)$parametersMAP)))

new_params=out_sample%>%
  group_by(variable_name)%>%
  summarise(value=median(value,na.rm=T))

new_params=merge(new_params,MAP_params)

new_params=left_join(by='variable_name',
                     data.frame(prior=prior$best,
                                lower=prior$lower,
                                upper=prior$upper,
                                variable_name=bayesianSetup$names
                     ),
                     data.frame(new_params))%>%
  mutate(post=ifelse(is.na(value),prior,value),
         value=NULL)

out_sample=merge(out_sample,new_params)

posterior_quantiles=out_sample%>%
  group_by(variable_name) %>% 
  summarise(post_lower=quantile(value,c(0.05)),
            post_upper=quantile(value,c(0.95)))

posterior_quantiles$prior_lower=prior$lower
posterior_quantiles$prior_upper=prior$upper

posterior_quantiles=merge(posterior_quantiles,
                          dplyr::select(new_params,variable_name,post,MAP),
                          by='variable_name')
posterior_quantiles=posterior_quantiles%>%
  mutate(uncertanty_reduction=100*(1-((post_upper-post_lower)/(prior_upper-prior_lower))),
         post_sd=((post-post_lower)+(post_upper-post))/2)

ggplot(posterior_quantiles,
              aes(x=variable_name,
                  y=post
              ))+
  facet_wrap(~variable_name, scales= "free")+
  geom_hline(aes(yintercept=prior_lower), linewidth = 0.1)+
  geom_hline(aes(yintercept=prior_upper), linewidth = 0.1)+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=post_lower,ymax=post_upper))


densityC <- ggplot(out_sample%>%
                     filter(
                       #variable_name=='lue_51'
                     ),
                   aes(x=value))+
  geom_density(aes(fill=variable_name),adjust=10,show.legend = F)+
  geom_vline(aes(xintercept=lower),size=0)+
  geom_vline(aes(xintercept=upper),size=0)+
  geom_vline(aes(xintercept=MAP),linetype='dashed')+
  geom_vline(aes(xintercept=post),linetype='dotted')+
  ggtitle("Ambient model")+
  facet_wrap(vars(variable_name),scales='free')
densityC

posterior_quantilesC=posterior_quantiles
final_parsC=posterior_quantiles$MAP


# Run MCMC Warming --------------------------------------------------------
prior=createUniformPrior(lower=c(0.002,0.0006,0.0002,
                                 10, 0.0025,0.0001),
                         upper=c(0.2,0.1,0.1,
                                 235,0.5,0.1))
settings = list(iterations = 500000 ,message=T)
bayesianSetup <- createBayesianSetup(get_likelihood, 
                                     prior=prior, 
                                     names = c('k1','k2','k3','input', 'a21','a23')
)

dat <-  data.frame(subset(datComb, year == 2009 | treatment == "Warming") %>% 
                     dplyr::select(-treatment) %>% 
                     mutate(C_stock = C_stock*1000,
                            C_stock_sd = C_stock_sd*1000))

# outW <- runMCMC(bayesianSetup = bayesianSetup,
#                 sampler = "DREAMzs",
#                 settings = settings)
# saveRDS(outW, file=paste0(MCMC_path,'mcmc_Out_Warming',format(Sys.Date(), "%Y%m%d%S"), Sys.time(),'.RDS'))

outW <- readRDS("~/Documents/NAU/SchuurLab/RCode/SoilR/MCMCOutput/mcmc_Out_Warming20250508002025-05-08 08:49:20.122385.RDS")
out <- outW

plot(out)
burn_in=150000
out_sample=reshape2::melt(variable.name='variable_name',
                          data.frame(getSample(out,start=burn_in)))
MAP_params=data.frame(MAP=MAP(out,start=burn_in)$parametersMAP,
                      variable_name=rownames(data.frame(MAP(out)$parametersMAP)))

new_params=out_sample%>%
  group_by(variable_name)%>%
  summarise(value=median(value,na.rm=T))

new_params=merge(new_params,MAP_params)

new_params=left_join(by='variable_name',
                     data.frame(prior=prior$best,
                                lower=prior$lower,
                                upper=prior$upper,
                                variable_name=bayesianSetup$names
                     ),
                     data.frame(new_params))%>%
  mutate(post=ifelse(is.na(value),prior,value),
         value=NULL)

out_sample=merge(out_sample,new_params)

posterior_quantiles=out_sample%>%
  group_by(variable_name) %>% 
  summarise(post_lower=quantile(value,c(0.05)),
            post_upper=quantile(value,c(0.95)))

posterior_quantiles$prior_lower=prior$lower
posterior_quantiles$prior_upper=prior$upper

posterior_quantiles=merge(posterior_quantiles,
                          dplyr::select(new_params,variable_name,post,MAP),
                          by='variable_name')
posterior_quantiles=posterior_quantiles%>%
  mutate(uncertanty_reduction=100*(1-((post_upper-post_lower)/(prior_upper-prior_lower))),
         post_sd=((post-post_lower)+(post_upper-post))/2)

pQW <- ggplot(posterior_quantiles,
              aes(x=variable_name,
                  y=post
              ))+
  facet_wrap(~variable_name, scales= "free")+
  geom_hline(aes(yintercept=prior_lower),size=0.1)+
  geom_hline(aes(yintercept=prior_upper),size=0.1)+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=post_lower,ymax=post_upper))

pQW

densityW <- ggplot(out_sample%>%
                     filter(
                       #variable_name=='lue_51'
                     ),
                   aes(x=value))+
  geom_density(aes(fill=variable_name),adjust=10,show.legend = F)+
  geom_vline(aes(xintercept=lower),size=0)+
  geom_vline(aes(xintercept=upper),size=0)+
  geom_vline(aes(xintercept=MAP),linetype='dashed')+
  geom_vline(aes(xintercept=post),linetype='dotted')+
  ggtitle("Warming model")+
  facet_wrap(vars(variable_name),scales='free')
densityW

posterior_quantilesW=posterior_quantiles
posterior_quantilesW$treatment = "Warming"
final_parsW=posterior_quantiles$MAP

# Plot final posterior draws ----------------------------------------------
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = alpha(c("#9AA582FF", "#657359FF", "#8B775FFF", "#6B472CFF"), alpha = 0.75)))


treatmentPosteriors <- posterior_quantilesC %>% 
  mutate(treatment = "Ambient") %>% 
  rbind(posterior_quantilesW) %>% 
  mutate(variable_name = factor(variable_name, levels = c("input","k1", "k2", "k3", "a21", "a23"),
                                labels = c(
                                  "C inputs (gC m<sup>−2</sup> yr<sup>−1</sup>)",
                                  "Surface organic k (yr<sup>−1</sup>)",
                                  "Deep organic k (yr<sup>−1</sup>)",
                                  "Mineral k (yr<sup>−1</sup>)",
                                  "Surf. to deep org. C transfer (yr<sup>−1</sup>)",
                                  "Deep org. to min. C transfer (yr<sup>−1</sup>)"
                                )))

decomp <- ggplot(subset(treatmentPosteriors, variable_name %in% c("Surface organic k (yr<sup>−1</sup>)",
                                                                  "Deep organic k (yr<sup>−1</sup>)",
                                                                  "Mineral k (yr<sup>−1</sup>)")),
                 aes(x=treatment,
                     y = MAP,
                     #y = post,
                     fill = treatment
                 ))+
  #facet_wrap(~variable_name, scales= "free", nrow = 3)+
  ggh4x::facet_wrap2(~variable_name, nrow = 3, scales = "free_y", strip = strip) +
  geom_bar(stat='identity', color = "black", position = position_dodge(), alpha = 0.75)+
  geom_errorbar(aes(x = treatment, ymin=post_lower,ymax=post_upper), width = .5)+
  scale_fill_manual(values = treatPal[-1])+
  ylab(element_blank())+xlab(element_blank())+
  theme(legend.position = "none", strip.text = ggtext::element_markdown())

cmove <- ggplot(subset(treatmentPosteriors, !(variable_name %in% c("Surface organic k (yr<sup>−1</sup>)",
                                                                   "Deep organic k (yr<sup>−1</sup>)",
                                                                   "Mineral k (yr<sup>−1</sup>)"))),
                aes(x=treatment,
                    y = MAP,
                    #y = post,
                    fill = treatment
                ))+
  facet_wrap(~variable_name, scales= "free_y", nrow = 3)+
  geom_bar(stat='identity', color = "black", position = position_dodge(), alpha = 0.75)+
  geom_errorbar(aes(x = treatment, ymin=post_lower,ymax=post_upper), width = .5)+
  scale_fill_manual(values = treatPal[-1])+
  ylab("Parameter estimate")+xlab(element_blank())+
  theme(legend.position = "none", strip.text = ggtext::element_markdown())

finPar <- cmove | decomp + patchwork::plot_layout(guides = "collect")
finPar

# Plot final model variables ----------------------------------------------
dat <- datComb %>% 
  mutate(C_stock = C_stock*1000,
         C_stock_sd = C_stock_sd*1000)

model_outputC=ThreepSeriesModel14(
  t=years,ks=c(k1=final_parsC[4], k2=final_parsC[5], k3=final_parsC[6]),
  C0=dat$C_stock[dat$year==2009], F0_Delta14C=dat$C14[dat$year==2009],
  In=final_parsC[3], a21=final_parsC[1], a32=final_parsC[2],
  inputFc=atmIn,pass=T
)

C14C=getF14(model_outputC)#[years%in%c(2009,2022),]
colnames(C14C) <- c("so.C14", "do.C14", "min.C14")
C_stockC=getC(model_outputC)#[years%in%c(2009,2022),]
colnames(C_stockC) <- c("so.CStock", "do.CStock", "min.CStock")

model_outputW=ThreepSeriesModel14(
  t=years,ks=c(k1=final_parsW[4], k2=final_parsW[5], k3=final_parsW[6]),
  C0=dat$C_stock[dat$year==2009], F0_Delta14C=dat$C14[dat$year==2009],
  In=final_parsW[3], a21=final_parsW[1], a32=final_parsW[2],
  inputFc=atmIn,pass=T
)

C14W=getF14(model_outputW)#[years%in%c(2009,2022),]
colnames(C14W) <- c("so.C14", "do.C14", "min.C14")
C_stockW=getC(model_outputW)#[years%in%c(2009,2022),]
colnames(C_stockW) <- c("so.CStock", "do.CStock", "min.CStock")

finalDat <- data.frame(cbind(C14C, C_stockC, years), treatment = "Control") %>% 
  rbind(data.frame(cbind(C14W, C_stockW, years), treatment = "Warming")) %>% 
  left_join(atmIn %>% rename(years = YEAR, atm.C14 = Atmosphere), by = "years") %>% 
  tidyr::pivot_longer(cols = -c(years, treatment),
                      names_to = c("layer", ".value"), 
                      names_sep = "\\.") %>% 
  mutate(layer = factor(layer, levels = c("so", "do", "min"), 
                        labels = c("Surface organic (~0-15cm)", "Deep organic (~15-35cm)", 
                                   "Mineral (~35-55cm)")))

datError <- subset(datComb, year == 2009 | treatment == "Control") %>% 
  mutate(treatment = "Control") %>% 
  rbind(mutate(subset(datComb, year == 2009 | treatment == "Warming"), treatment = "Warming")) %>% 
  mutate(layer = factor(layer, levels = c("so", "do", "min"), 
                        labels = c("Surface organic (~0-15cm)", "Deep organic (~15-35cm)", 
                                   "Mineral (~35-55cm)")))
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = alpha(c("#9AA582FF", "#657359FF", "#8B775FFF", "#6B472CFF"), alpha = 0.75)))

cStockByLayer <- ggplot()+
  geom_line(data = subset(finalDat, layer != "atm"), 
            aes(x = years, y = CStock/1000, group = treatment, color = treatment), 
            size = 1.5, alpha = 0.75) +
  geom_point(data = datError, aes(x = year, y = (C_stock), color = treatment), size = 4, show.legend = FALSE)+
  geom_errorbar(data = datError, aes(x = year, ymin = (C_stock - C_stock_sd), ymax = (C_stock + C_stock_sd), color = treatment), 
                size = 1, width = 0.5, alpha = 0.75, show.legend = FALSE)+
  scale_color_manual(values = treatPal[-1], name = "Treatment")+
  ylab(expression(paste("C stock ", "(kg m"^-2, ")")))+xlab("Years")+
  theme(legend.position = "none")+
  ggh4x::facet_wrap2(~layer, nrow = 3, scales = "free_y", strip = strip) 
cStockByLayer

c14ByLayer <- ggplot()+
  geom_line(data = subset(finalDat, !is.na(C14) & !is.na(layer)), aes(x = years, y = C14, group = treatment, color = treatment), size = 1.5, alpha = 0.75)+
  #geom_line(data = atmCip, aes(x = YEAR, y = Atmosphere, color = "atm"), size = 1.5, alpha = 0.75, color = "#808080")+
  geom_point(data = datError, aes(x = year, y = C14, color = treatment), size = 4)+
  geom_errorbar(data = datError, aes(x = year, ymin = C14 - C_14_sd, ymax = C14 + C_14_sd, color = treatment), size = 1, width = 0.5, alpha = 0.75)+  
  scale_color_manual(values = treatPal[-1], name = "Treatment")+
  ylab(expression(paste(Delta^14,"C (‰)")))+ xlab("Years")+
  theme(legend.position = "none")+
  ggh4x::facet_wrap2(~layer, nrow = 3, scales = "free_y", strip = strip) 
c14ByLayer

fin <- ((c14ByLayer | cStockByLayer))
fin 

finPar | fin

#final 14C modeling figure

((c14ByLayer | cStockByLayer) | cmove | decomp) + 
  patchwork::plot_layout(guides = "collect", widths = c(3.2,3.2, 2.95, 2.75)) + 
  patchwork::plot_annotation(tag_levels = list(c('a)', 'b)','c)','' ), '2')) & theme(plot.tag = element_text(size = 20))
#ggsave(paste0(manDirFin, "mcmcOutput.png"), width = 40, height = 25, units = "cm", dpi = 300)



# Supplementary MCMC figures ----------------------------------------------

densityC / densityW
#ggsave(paste0(manDirFin, "mcmcSupp.png"), width = 30, height = 20, units = "cm", dpi = 300)
