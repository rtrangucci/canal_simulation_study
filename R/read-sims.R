source("R/summarise_simulations.R")
library(ggplot2)

all_res <- read_files("mez_sim_res")
names(all_res)

 all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_x=mean(mse_rho,na.rm=T),sd_mse = sd(mse_rho,na.rm=T)/sqrt(50), bias = mean(rho,na.rm=T) - (0.1),
            sd = mean(sd_rho^2,na.rm=T),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform'))
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat) %>% xtable::xtable(digits = 5)

rho_sum <- all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_x=mean(bias_rho^2,na.rm=T),
            sd_mse = sd(bias_rho^2,na.rm=T)/sqrt(n()),
            bias = mean(bias_rho,na.rm=T),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform')),
    upper = mse_x + 2*sd_mse,
    lower = mse_x - 2*sd_mse
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat)

filter(rho_sum, J > 100) %>% 
  ggplot(aes(x = J, y = mse_x, colour = interaction(I,Distribution), group = interaction(I,Distribution))) +
  geom_line() +
  xlab('Number of observations') + ylab(bquote('MSE '*rho)) +
  scale_color_discrete(name = "Spatial distribution") + theme_bw()

all_res %>%
  filter(n > 100) |>
  group_by(G, n, spat) %>%
  summarise(
    mse_x = mean(bias_alpha^2, na.rm = T),
    sd_mse = sd(bias_alpha^2, na.rm = T) / sqrt(n()),
    bias = mean(bias_alpha, na.rm = T),
    tot_n = tot_n[1]
  ) %>%
  mutate(
    spat = plyr::mapvalues(spat, c("c", "u"), c("Clustered", "Uniform")),
    upper = mse_x + 2 * sd_mse,
    lower = mse_x - 2 * sd_mse
  ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat) |>
  ggplot(aes(x = J, y = mse_x, colour = interaction(I, Distribution), group = interaction(I, Distribution))) +
  geom_line() +
  xlab("Number of observations") +
  ylab(bquote("MSE " * alpha)) +
  scale_color_discrete(name = "Spatial distribution") +
  theme_bw()

filter(rho_sum, I == 100) %>%
  ggplot(aes(x = N, y = mse_x, group = Distribution)) +
  geom_line(aes(colour = Distribution)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 100) +
  xlab("Number of observations") +
  ylab(bquote("MSE " * rho)) +
  scale_color_discrete(name = "Spatial distribution") +
  theme_bw()

all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_x=mean(mse_beta,na.rm=T),sd_mse = sd(mse_beta,na.rm=T)/sqrt(50), bias = mean(beta,na.rm=T) - (-0.15),
            sd = mean(sd_alpha^2,na.rm=T),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform'))
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat, MSE = mse_x, `SD MSE` = sd_mse, Bias = bias, Var = sd) %>% xtable(digits = 5)

imse_plot <- all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_x=mean(classical_imse_x,na.rm=T),sd_mse_x = sd(classical_imse_x,na.rm=T)/sqrt(n()),
            mse_x2=mean(classical_imse_x2,na.rm=T),sd_mse_x2 = sd(classical_imse_x2,na.rm=T)/sqrt(n()),
            mse_y=mean(classical_imse_x2,na.rm=T),sd_mse_y = sd(classical_imse_y,na.rm=T)/sqrt(n()),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform'))
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat,
         `$MSE_x$` = mse_x, `SD $MSE_x$` = sd_mse_x,
         `$MSE_{x_2}$` = mse_x2, `SD $MSE_{x_2}$` = sd_mse_x2,
         `$MSE_y$` = mse_y, `SD $MSE_y$` = sd_mse_y,
         ) %>% xtable::xtable(digits = 3)

imse_plot <- all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_x=mean(classical_imse_x,na.rm=T),sd_mse_x = sd(classical_imse_x,na.rm=T)/sqrt(n()),
            mse_x2=mean(classical_imse_x2,na.rm=T),sd_mse_x2 = sd(classical_imse_x2,na.rm=T)/sqrt(n()),
            mse_y=mean(classical_imse_x2,na.rm=T),sd_mse_y = sd(classical_imse_y,na.rm=T)/sqrt(n()),
            var_x=mean(classical_ivar_x), sd_var_x = sd(classical_ivar_x)/sqrt(n()),
            var_y=mean(classical_ivar_y), sd_var_y = sd(classical_ivar_y)/sqrt(n()),
            var_x2=mean(classical_ivar_x2), sd_var_x2 = sd(classical_ivar_x2)/sqrt(n()),
            bias_x=mean(classical_ibias_x), sd_bias_x = sd(classical_ibias_x)/sqrt(n()),
            bias_y=mean(classical_ibias_y), sd_bias_y = sd(classical_ibias_y)/sqrt(n()),
            bias_x2=mean(classical_ibias_x2), sd_bias_x2 = sd(classical_ibias_x2)/sqrt(n()),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform')),
    upper_mse_x = mse_x + 2*sd_mse_x,
    lower_mse_x = mse_x - 2*sd_mse_x,
    upper_mse_x2 = mse_x2 + 2*sd_mse_x2,
    lower_mse_x2 = mse_x2 - 2*sd_mse_x2,
    upper_mse_y = mse_y + 2*sd_mse_y,
    lower_mse_y = mse_y - 2*sd_mse_y,
    upper_var_x = var_x + 2*sd_var_x,
    lower_var_x = var_x - 2*sd_var_x,
    upper_var_x2 = var_x2 + 2*sd_var_x2,
    lower_var_x2 = var_x2 - 2*sd_var_x2,
    upper_var_y = var_y + 2*sd_var_y,
    lower_var_y = var_y - 2*sd_var_y,
    upper_bias_x = bias_x + 2*sd_bias_x,
    lower_bias_x = bias_x - 2*sd_bias_x,
    upper_bias_x2 = bias_x2 + 2*sd_bias_x2,
    lower_bias_x2 = bias_x2 - 2*sd_bias_x2,
    upper_bias_y = bias_y + 2*sd_bias_y,
    lower_bias_y = bias_y - 2*sd_bias_y
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat)

mse_plot <- all_res %>%
  group_by(G,n,spat) %>%
  summarise(mse_rho=mean(bias_rho^2,na.rm=T),
            mse_alpha = mean(bias_alpha^2, na.rm = T),
            mse_beta = mean(bias_beta^2, na.rm = T),
            se_mse_rho=sd(bias_rho^2,na.rm=T)/sqrt(n()),
            se_mse_alpha = sd(bias_alpha^2, na.rm = T)/sqrt(n()),
            se_mse_beta = sd(bias_beta^2, na.rm = T)/sqrt(n()),
            bias_rho = mean(bias_rho,na.rm=T), bias_alpha = mean(bias_alpha, na.rm=T),
            bias_beta = mean(bias_beta, na.rm=T),
            sd_bias_rho = sd(bias_rho,na.rm=T)/sqrt(n()), sd_bias_alpha = sd(bias_alpha,na.rm=T)/sqrt(n()),
            sd_bias_beta = sd(bias_beta,na.rm=T)/sqrt(n()),
            var_rho = mean(sd_rho^2,na.rm=T),
            var_alpha = mean(sd_alpha^2,na.rm=T),
            var_beta = mean(sd_beta^2,na.rm=T),
            se_var_rho = sd(sd_rho^2,na.rm=T)/sqrt(n()),
            se_var_alpha = sd(sd_alpha^2,na.rm=T)/sqrt(n()),
            se_var_beta = sd(sd_beta^2,na.rm=T)/sqrt(n()),sd_rho = sd(mse_rho,na.rm=T)/sqrt(n()),
            se_alpha=mean(mse_alpha,na.rm=T),sd_alpha = sd(mse_alpha,na.rm=T)/sqrt(n()),
            se_beta=mean(mse_beta,na.rm=T),sd_beta = sd(mse_beta,na.rm=T)/sqrt(n()),
            cover_theta_50 = mean(pred_cover_50,na.rm=T),sd_cover_theta_50 = sd(pred_cover_50,na.rm=T)/sqrt(n()),
            cover_theta_80 = mean(pred_cover_80,na.rm=T),sd_cover_theta_80 = sd(pred_cover_80,na.rm=T)/sqrt(n()),
            cover_risk_50 = mean(risk_cover_50,na.rm=T),sd_cover_risk_50 = sd(risk_cover_50,na.rm=T)/sqrt(n()),
            cover_risk_80 = mean(risk_cover_80,na.rm=T),sd_cover_risk_80 = sd(risk_cover_80,na.rm=T)/sqrt(n()),
            bias_theta = mean(pred_bias,na.rm=T),sd_bias = sd(pred_bias,na.rm=T)/sqrt(n()),
            bias_risk = mean(risk_bias,na.rm=T),sd_risk_bias = sd(risk_bias,na.rm=T)/sqrt(n()),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform')),
    upper_alpha = se_alpha + 2*sd_alpha,
    lower_alpha = se_alpha - 2*sd_alpha,
    upper_rho = mse_rho + 2*se_mse_rho,
    lower_rho = mse_rho - 2*se_mse_rho,
    upper_beta = mse_beta + 2*se_mse_beta,
    lower_beta = mse_beta - 2*se_mse_beta,
    upper_cover_theta_50 = cover_theta_50 + 2*sd_cover_theta_50,
    lower_cover_theta_50 = cover_theta_50 - 2*sd_cover_theta_50,
    upper_cover_theta_80 = cover_theta_80 + 2*sd_cover_theta_80,
    lower_cover_theta_80 = cover_theta_80 - 2*sd_cover_theta_80,
    upper_cover_risk_50 = cover_risk_50 + 2*sd_cover_risk_50,
    lower_cover_risk_50 = cover_risk_50 - 2*sd_cover_risk_50,
    upper_cover_risk_80 = cover_risk_80 + 2*sd_cover_risk_80,
    lower_cover_risk_80 = cover_risk_80 - 2*sd_cover_risk_80,
    upper_bias = bias_theta + 2*sd_bias,
    lower_bias = bias_theta - 2*sd_bias,
    upper_risk_bias = bias_risk + 2*sd_risk_bias,
    lower_risk_bias = bias_risk - 2*sd_risk_bias,
    upper_bias_rho = bias_rho + 2*sd_bias_rho,
    lower_bias_rho = bias_rho - 2*sd_bias_rho,
    upper_bias_alpha = bias_alpha + 2*sd_bias_alpha,
    lower_bias_alpha = bias_alpha - 2*sd_bias_alpha,
    upper_bias_beta = bias_beta + 2*sd_bias_beta,
    lower_bias_beta = bias_beta - 2*sd_bias_beta,
    upper_var_rho = var_rho + 2*se_var_rho,
    lower_var_rho = var_rho - 2*se_var_rho,
    upper_var_alpha = var_alpha + 2*se_var_alpha,
    lower_var_alpha = var_alpha - 2*se_var_alpha,
    upper_var_beta = var_beta + 2*se_var_beta,
    lower_var_beta = var_beta - 2*se_var_beta,
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat)

rat_var_res <- all_res %>%
  group_by(G,n) %>% 
  summarise(
    mean_rho_c = mean(sd_rho[spat == 'c']^2,na.rm=T),
    mean_rho_u = mean(sd_rho[spat == 'u']^2,na.rm=T),
    cov_rho = cov(sd_rho[spat == 'c']^2,sd_rho[spat == 'u']^2, use = 'pairwise.complete'),
    rat_rho = mean_rho_c / mean_rho_u,
    sd_alt_rat = 1/sqrt(n())/mean_rho_u * 
      sqrt(rat_rho^2 * var(sd_rho[spat == 'u']^2,na.rm=T) + var(sd_rho[spat == 'c']^2,na.rm=T)),
    sd_rat_rho = rat_rho / sqrt(n())* 
      sqrt(var(sd_rho[spat == 'c']^2,na.rm=T)/mean_rho_c^2 + var(sd_rho[spat == 'u']^2,na.rm=T)/mean_rho_u^2),
    mean_alpha_c = mean(sd_alpha[spat == 'c']^2,na.rm=T),
    mean_alpha_u = mean(sd_alpha[spat == 'u']^2,na.rm=T),
    rat_alpha = mean_alpha_c / mean_alpha_u,
    sd_rat_alpha = rat_alpha * sqrt(var(sd_alpha[spat == 'c']^2,na.rm=T)/mean_alpha_c^2 + var(sd_alpha[spat == 'u']^2,na.rm=T)/mean_alpha_u^2) / sqrt(n()),
  ) %>% mutate(
    N = n * G,
    upper_rho = rat_rho + 2 * sd_rat_rho,
    lower_rho = rat_rho - 2 * sd_rat_rho,
    upper_alpha = rat_alpha + 2 * sd_rat_alpha,
    lower_alpha = rat_alpha - 2 * sd_rat_alpha,
  ) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

side_by_side <- function(g1, g2) {
  legend <- g_legend(g1)
  g3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(g1 + theme(legend.position="none"),
                           g2 + theme(legend.position="none"),
                           nrow=1),
               legend, nrow=2,heights=c(10, 1)) 

return(g3)
}

  ggplot(aes(x = J, y = mse_x, colour = interaction(I,Distribution), group = interaction(I,Distribution))) +


p1 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = J, y = mse_x,
             group = interaction(I,Distribution),
         colour = interaction(I,Distribution))) +
  geom_line() +
  xlab('Number of observations') + ylab(bquote('IMSE '*Lambda[x])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean IMSE '*Lambda[x])) + theme_bw() +
  guides(colour = guide_legend(nrow=1)) 


p2 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = mse_y, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_mse_y, ymax = upper_mse_y),width = 100) + 
  xlab('Number of observations') + ylab(bquote('IMSE '*Lambda[y])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean IMSE '*Lambda[y])) + theme_bw() +
  guides(colour = guide_legend(nrow=1))

p3 <- side_by_side(p1,p2)

p1 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_x, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_var_x, ymax = upper_var_x),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ivar '*Lambda[x])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ivar '*Lambda[x])) + theme_bw() +
  guides(colour = guide_legend(nrow=1)) 


p2 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_y, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_var_y, ymax = upper_var_y),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ivar '*Lambda[y])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ivar '*Lambda[y])) + theme_bw() +
  guides(colour = guide_legend(nrow=1))
p3 <- side_by_side(p1,p2)

p1 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_x, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias_x, ymax = upper_bias_x),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ibias '*Lambda[x])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ibias '*Lambda[x])) + theme_bw() +
  guides(colour = guide_legend(nrow=1)) 


p2 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_y, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias_y, ymax = upper_bias_y),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ibias '*Lambda[y])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ibias '*Lambda[y])) + theme_bw() +
  guides(colour = guide_legend(nrow=1))
p3 <- side_by_side(p1,p2)

p1 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_x2, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias_x2, ymax = upper_bias_x2),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ibias '*Lambda[x])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ibias '*Lambda[x])) + theme_bw() +
  guides(colour = guide_legend(nrow=1)) 


p2 <- filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_y, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias_y, ymax = upper_bias_y),width = 100) + 
  xlab('Number of observations') + ylab(bquote('Ibias '*Lambda[y])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean Ibias '*Lambda[y])) + theme_bw() +
  guides(colour = guide_legend(nrow=1))
p3 <- side_by_side(p1,p2)

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_rho, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  geom_errorbar(aes(ymin = lower_rho, ymax = upper_rho),width = 100, position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('MSE '*rho)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*rho)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_alpha, ymax = upper_alpha),width = 100) + 
  xlab('Number of observations') + ylab(bquote('MSE '*lambda)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_rho, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*rho)) +
  geom_errorbar(aes(ymin = lower_bias_rho, ymax = upper_bias_rho),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*rho)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*lambda)) +
  geom_errorbar(aes(ymin = lower_bias_alpha, ymax = upper_bias_alpha),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))  

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_rho, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*rho)) +
  geom_errorbar(aes(ymin = lower_var_rho, ymax = upper_var_rho),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*rho)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*lambda)) +
  geom_errorbar(aes(ymin = lower_var_alpha, ymax = upper_var_alpha),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p1 <- p1 + ylim(range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))
p2 <- p2 +  ylim(range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))  

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_rho, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*rho)) +
  geom_errorbar(aes(ymin = lower_var_rho, ymax = upper_var_rho),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*rho)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

with(rat_var_res %>% filter(n > 150),plot(N, rat_rho,type = 'l',ylim = c(0.8,2)))
with(rat_var_res %>% filter(n > 150),lines(N, rat_alpha,col='red'))

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*lambda)) +
  geom_errorbar(aes(ymin = lower_var_alpha, ymax = upper_var_alpha),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p1 <- p1 + ylim(range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))
p2 <- p2 +  ylim(range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))  

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_beta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*beta)) +
  geom_errorbar(aes(ymin = lower_var_beta, ymax = upper_var_beta),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*beta)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = var_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Var. '*lambda)) +
  geom_errorbar(aes(ymin = lower_var_alpha, ymax = upper_var_alpha),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean var '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                          p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))  

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_beta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*gamma)) +
  geom_errorbar(aes(ymin = lower_bias_beta, ymax = upper_bias_beta),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*gamma)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

plt_df <- reshape2::melt(rat_var_res,measure.vars = c('rat_rho','rat_alpha'))  
plt_df$upper <- c(plt_df$upper_rho[plt_df$variable == 'rat_rho'],
                  plt_df$upper_alpha[plt_df$variable == 'rat_alpha'])
plt_df$lower <- c(plt_df$lower_rho[plt_df$variable == 'rat_rho'],
                  plt_df$lower_alpha[plt_df$variable == 'rat_alpha'])

plt_df %>%
  filter(n > 150) %>%
  ggplot(aes(x = N, y = value,
             group=variable)) + 
  geom_line(aes(colour=variable)) + 
  geom_errorbar(aes(ymax=upper,ymin = lower),width=100) +
  theme_bw() + ylab('Variance ratio') + xlab('Number of observations') + scale_color_discrete(name = 'Parameter', labels = c(as.expression(bquote(rho)),as.expression(bquote(lambda)))) -> p_rat

ggsave(p_rat, filename = 'prelim/test-rat-var-rho-lambda.pdf', width = 8.51, height = 4.94,units = 'in')

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

all_res %>%
  filter(n > 150,spat == 'c') %>%
  mutate(N = n * G) %>%
  ggplot(aes(y = sd_rho^2, x = sd_alpha^2, group = as.factor(N))) +
  geom_point(aes(colour=as.factor(N))) +
  theme_bw() + xlab(bquote('Posterior variance '*lambda)) + ylab(bquote('Posterior variance '*rho)) +
  scale_colour_manual('Number of observations', values = cbp1) + guides(colour = guide_legend(nrow=1)) +
  ggtitle('Clustered households') -> p_c_var

all_res %>%
  filter(n > 150,spat == 'u') %>%
  mutate(N = n * G) %>%
  ggplot(aes(y = sd_rho^2, x = sd_alpha^2, group = as.factor(N))) +
  geom_point(aes(colour=as.factor(N))) +
  theme_bw() + xlab(bquote('Posterior variance '*lambda)) + ylab(bquote('Posterior variance '*rho)) +
  scale_colour_manual('Number of observations', values = cbp1) + guides(colour = guide_legend(nrow=1)) + 
  ggtitle('Uniformly distributed households') -> p_u_var

get_axis_range <- function(g1, g2, axis = 'x') {
  if (axis == 'x')
    return(range(c(ggplot_build(g1)$layout$panel_scales_x[[1]]$range$range,ggplot_build(g2)$layout$panel_scales_x[[1]]$range$range)))
  else
    return(range(c(ggplot_build(g1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(g2)$layout$panel_scales_y[[1]]$range$range)))
}


y_range <- get_axis_range(p_c_var, p_u_var, 'y')
x_range <- get_axis_range(p_c_var, p_u_var, 'x')
p_c_var <- p_c_var + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = y_range) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), limits = x_range)  
p_u_var <- p_u_var + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = y_range) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE),limits = x_range)  
p_comb <- side_by_side(p_c_var, p_u_var)
ggsave(p_comb, file = 'prelim/var-rho-lambda-scatter.pdf', width = 8.51, height = 4.94, units = 'in')

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*lambda)) +
  geom_errorbar(aes(ymin = lower_bias_alpha, ymax = upper_bias_alpha),width = 100, position=position_dodge(width = 200)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))  
p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_beta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_beta, ymax = upper_beta),width = 100, position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('MSE '*gamma)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*gamma)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_alpha, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_alpha, ymax = upper_alpha),width = 100) + 
  xlab('Number of observations') + ylab(bquote('MSE '*lambda)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*lambda)) +
  guides(colour = guide_legend(nrow=1)) + theme_bw()
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_rho, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_rho, ymax = upper_rho),width = 100, position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('MSE '*rho)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*rho))

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = se_beta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_beta, ymax = upper_beta),width = 100, position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('MSE '*gamma)) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('MSE '*gamma))
gridExtra::grid.arrange(p1,p2, nrow = 1, ncol = 2)

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = cover_theta_50, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_cover_theta_50, ymax = upper_cover_theta_50),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Coverage '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean coverage '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0.5,colour = 'black',linetype = 'dotted')

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_theta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias, ymax = upper_bias),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0,colour = 'black',linetype = 'dotted')
p3 <- side_by_side(p1,p2)

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = cover_theta_80, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_cover_theta_80, ymax = upper_cover_theta_80),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Coverage '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean coverage '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0.8,colour = 'black',linetype = 'dotted')

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_theta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias, ymax = upper_bias),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0,colour = 'black',linetype = 'dotted')
p3 <- side_by_side(p1,p2)

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_risk, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_risk_bias, ymax = upper_risk_bias),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0,colour = 'black',linetype = 'dotted')

p1 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = cover_theta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_cover, ymax = upper_cover),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Coverage '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean coverage '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0.5,colour = 'black',linetype = 'dotted')

p2 <- filter(mse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = bias_theta, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_bias, ymax = upper_bias),width = 100,
                                                        position=position_dodge(width = 200)) + 
  xlab('Number of observations') + ylab(bquote('Bias '*theta[i])) +
  scale_color_discrete(name = "Spatial distribution") + ggtitle(bquote('Mean bias '*theta[i])) +
  guides(colour = guide_legend(nrow=1)) + theme_bw() + geom_hline(yintercept = 0,colour = 'black',linetype = 'dotted')
mylegend<-g_legend(p1)
p3 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))

filter(imse_plot, J > 150) %>% 
  ggplot(aes(x = N, y = mse_x2, group = Distribution)) +
  geom_line(aes(colour = Distribution)) + geom_errorbar(aes(ymin = lower_x2, ymax = upper_x2),width = 100) + 
  xlab('Number of observations') + ylab(bquote('MSE '*Lambda[x[2]])) +
  scale_color_discrete(name = "Spatial distribution")

all_res %>%
  group_by(G,n,spat) %>%
  summarise(cover=mean(pred_cover,na.rm=T),sd_cover = sd(pred_cover,na.rm=T)/sqrt(50),
            bias=mean(pred_bias,na.rm=T),sd_bias = sd(pred_bias,na.rm=T)/sqrt(50),
            tot_n = tot_n[1]) %>% 
  mutate(
    spat = plyr::mapvalues(spat, c('c','u'), c('Clustered','Uniform'))
    ) %>%
  rename(I = G, J = n, N = tot_n, Distribution = spat,
         Cover = cover, `SD Cover` = sd_cover,
         Bias = bias, `SD Bias` = sd_bias
         ) %>% xtable(digits = 5) %>% print(include.rownames = FALSE) 
