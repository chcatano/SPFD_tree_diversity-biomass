## Analysis script for structural equation models 
## β-diversity as a driver of forest biomass across spatial scales
## Jacqueline C. Reu, Christopher P. Catano, Marko J. Spasojevic, Jonathan A. Myers


rm(list = ls())

# Check for and install/load required packages
for (package in c('car', 'tidyverse', 'piecewiseSEM', 'lattice', 'ggpubr', 'viridis')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. IMPORT DATA ####
my.data <- read.csv("data/trcp_bef_data_final.csv")
my.data$scale <- as.character(my.data$scale)
# log-transform biomass to linearize relationships
my.data$ln_biom <- log(my.data$mean_biomass)

#### 2. Fit Piecewise SEMs ####
# beta diversity path coefficients
beta.coefs <- my.data %>% group_by(scale) %>%
  do(coefs(psem(lm(ln_biom ~ mean_betaPIE + mean_env, data = .), 
                    lm(mean_betaPIE ~ mean_env, data = .))))

# gamma diversity path coefficients
gamma.coefs <- my.data %>% group_by(scale) %>%
  do(coefs(psem(lm(ln_biom ~ mean_gPIE + mean_env, data = .), 
                lm(mean_gPIE ~ mean_env, data = .))))

# alpha diversity path coefficients
alpha.coefs <- my.data %>% group_by(scale) %>%
  do(coefs(psem(lm(ln_biom ~ mean_aPIE + mean_env, data = .), 
                lm(mean_aPIE ~ mean_env, data = .))))


#### 3. Plot results (Fig. 3 in main manuscript) ####
# Set figure aesthetics
theme_set(theme_bw() +  
            theme(legend.position = "bottom", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), 
                  plot.title = element_text(hjust = 0.5, size = 12), 
                  text = element_text(size = 12), axis.text = element_text(size = 10)))

# beta diversity path coefficients
beta.coefs$scale <- as.numeric(beta.coefs$scale)
beta.coefs$scale <- beta.coefs$scale/100
final.beta <- beta.coefs %>%
  unite('path', Response:Predictor, remove = F)

(scale.beta <- ggplot(final.beta[final.beta$Response != "mean_betaPIE", ], aes(x = scale, y = Std.Estimate, 
                                          group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Std. path coefficient") +
    ylim(c(-0.1,0.8))+
    ggtitle("Beta-diversity") +
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
        scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#440154FF", "#FDE725FF", "#238A8DFF")) +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
)


# alpha diversity path coefficients
alpha.coefs$scale <- as.numeric(alpha.coefs$scale)
alpha.coefs$scale <- alpha.coefs$scale/100
final.alpha <- alpha.coefs %>%
  unite('path', Response:Predictor, remove = F)

(scale.alpha <- ggplot(final.alpha[final.alpha$Response != "mean_aPIE", ], aes(x = scale, y = Std.Estimate, 
                                      group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Std. path coefficient") +
    ylim(c(-0.1,0.8))+
    ggtitle("Alpha-diversity") +
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
    scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#440154FF", "#FDE725FF", "#238A8DFF")) +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
)


# gamma diversity path coefficients
gamma.coefs$scale <- as.numeric(gamma.coefs$scale)
gamma.coefs$scale <- gamma.coefs$scale/100
final.gamma <- gamma.coefs %>%
  unite('path', Response:Predictor, remove = F)

(scale.gamma <- ggplot(final.gamma[final.gamma$Response != "mean_gPIE", ], aes(x = scale, y = Std.Estimate, 
                                      group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Biomass response \n (Std. path coefficient)") +
    ylim(c(-0.1,0.8))+
    ggtitle("Gamma-diversity") +
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
    scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#FDE725FF", "#440154FF", "#238A8DFF")) +
    theme(legend.position = "none", 
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
)


# effects of env on gamma
(scale.gamma.ind <- ggplot(final.gamma[final.gamma$Response == "mean_gPIE", ], aes(x = scale, y = Std.Estimate, 
                                                                                   group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Diversity response \n (Std. path coefficient)") +
    ylim(c(-0.15,0.8))+
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
    scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#238A8DFF")) +
    theme(legend.position = "none", axis.title.x=element_blank(),
          axis.title.y=element_blank())
)


# effects of env on beta
(scale.beta.ind <- ggplot(final.beta[final.beta$Response == "mean_betaPIE", ], aes(x = scale, y = Std.Estimate, 
                                                                                   group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Std. path coefficient") +
    ylim(c(-0.15,0.8))+
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
    scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#238A8DFF")) +
    theme(legend.position = "none", axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
)


# effects of env on alpha
(scale.alpha.ind <- ggplot(final.alpha[final.alpha$Response == "mean_aPIE", ], aes(x = scale, y = Std.Estimate, 
                                                                                   group = path, color = path)) + 
    geom_point(size = 2, alpha = 0.6) + 
    geom_smooth(method = "loess", se = F) +
    xlab(expression(paste("Area (m"^2, " x 100)"))) +
    ylab("Std. path coefficient") +
    ylim(c(-0.15,0.8))+
    geom_hline(yintercept = 0, col = "dark gray", linetype = "dashed") +
    scale_x_continuous(breaks = c(4, 16, 36, 64, 100, 144)) +
    scale_color_manual(values = c("#238A8DFF")) +
    theme(legend.position = "none", axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
)


# output figure for SEM paths vs scale 
#jpeg(file = "....jpeg", width = 8, height = 6, units = 'in', res = 1000)
fig <- ggarrange(scale.gamma, scale.beta, scale.alpha, 
  scale.gamma.ind, scale.beta.ind, scale.alpha.ind,
         ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom", align = "hv")
annotate_figure(fig,
                bottom = text_grob(expression(paste("Area (m"^2, " x 100)"))))
#dev.off()