# Numerical considerations for ICS

# Libraries ----
## Data -----
library(MASS) # contains crabs dataset
library(VGAM) #rlaplace simulation
load("data/HTP2.rda")
load("data/HTP3.rda")

## Processing -----
library(magrittr) # pipe
library(janitor) # clean_names
library(data.table) #rbindlist

## Methods ----
library(ICS)
library(moments)

## Plot ----
library(ggplot2)
library(GGally) ## ggpairs
library(ggthemes) ## scale_colour_colorblind

## Customized ICS_QR functions
source("R/ICS_utilities.R")

## Initialization -----
alpha_all =  c(-1, -0.5, 0.5, 1)
set.seed(20212022)



# Applications of ICS for clustering -----

## Mixture of two normal distributions ----
### Simulation of data ----
p = 4
Ntot = 10000
eps = 0.10
delta = 6

# Simulation of mixture of two Gaussian distributions
X_ini <- simulations(n = Ntot*(1-eps), n_outlier = Ntot*eps, p = p,
                     sigma = diag(1,p), mean =  rep(1,p),
                     sigma_outlier = diag(1,p),
                     mean_outlier =  c(delta,rep(1,p-1)))

# Simulations of vectors c_k
k_all <- 0:30
mean_all_save <- sapply(k_all, function(k){
  k_sim(k = k, p=p)
})
colnames(mean_all_save) <- paste0("K", k_all)

# Description
# Initial data
png("figures/simulations_ggpairs.png", width=900, height=900)
GGally::ggpairs(data.frame(X_ini, group = c(rep("Group 0",Ntot*(1-eps)),rep("Group 1", Ntot*eps))), aes(colour = group,
                        alpha = 0.4),
                columns = 1:p,
                upper = list(continuous = wrap("cor", size = 5)),
                title = "Mixture of two Gaussian distributions") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1),
                     guide = guide_axis(angle = 45) ) +

  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
dev.off()




# Scaled data
k = 8
X <-  sweep(X_ini, 2, mean_all_save[,k+1], "*")
print(colMeans(X))
print(format(kappa(X, exact = TRUE), scientific = TRUE, digits = 2))

png(paste0("figures/simulations_ggpairs_k", k,".png"),width=900, height=900)
GGally::ggpairs(data.frame(X,
                                group = c(rep("Group 0",Ntot*(1-eps)),rep("Group 1", Ntot*eps))),
                     aes(colour = group, alpha = 0.4),
                upper = list(continuous = wrap("cor", size = 5)),
                columns = 1:p,
                title = bquote("Mixture of two Gaussian distributions"~kappa~"~"~ 10^ .(k)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1),
                     guide = guide_axis(angle = 45)) +

  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20, lineheight = 0.3))
dev.off()



### Computation of ICS eigenvalues by different methods ----
all_methods <- c("ICSEigen", "ICSQR", "ICSTheory")

lapply(alpha_all, function(alpha){

  eigenvalues_all <- sapply(k_all+1, function(k){
    eigenvalues_all_methods(X_ini, k, alpha, cf = 1)
  })
  eigenvalues_all_dt <- data.frame(eigenvalues_all)

  dt_eigenvalues_plot(all_methods, eigenvalues_all_dt, nb = ncol(X_ini),
                      alpha, df_name = "simulations_eigenvalues",
                      save = TRUE, k_all = k_all)
})


## Crabs data set ----
### Preprocessing ----
data(crabs)

# LOG transformation
X = apply(as.matrix(crabs[,-c(1:3)]), 2, log)
dim(X)
rownames(X) <- 1:nrow(X)

# format(kappa(X, exact = TRUE), scientific = TRUE, digits = 2)

df <- data.frame(sp = crabs$sp,
                 sex = crabs$sex,
                 group = paste(crabs$sp, crabs$sex, sep = "_"),
                 X)
### Initial data ----
pdf("figures/crabs_log_ggpairs.pdf")
GGally::ggpairs(df, aes(colour = group,
                        alpha = 0.4),
                columns = 4:ncol(crabs),
                title = "Log crabs data")
dev.off()

### ICS - QR
sapply(alpha_all, function(alpha){
  res_ICS_QR <- ICS_QR_full_rank(Xt = X, alpha = alpha, cf = 1)
  save_ggpairs(crabs, res_ICS_QR, df_name = "crabs_log", save = TRUE)
})



### Scaled data ----
# scale (k=6 OK but not k=7)
mean_all <- k_sim(k = 6 , p = ncol(X))
X = apply(as.matrix(crabs[,-c(1:3)]), 2, log)
X <-  sweep(X, 2, mean_all, "*")
print(colMeans(X))
print(format(kappa(X, exact = TRUE), scientific = TRUE,
             digits = 2))

# ICS - ICS_EIGEN
res_ICS_EIGEN <- ICS::ics2(X, S1 = MeanCov,
                           S2 = MeanCovW, S2args = list(alpha = alpha, cf = 1))



# ICS - QR
sapply(alpha_all, function(alpha){
  res_ICS_QR <- ICS_QR_full_rank(Xt = X, alpha = alpha, cf = 1)
  save_ggpairs(crabs, res_ICS_QR, df_name = "crabs_log_scale", save = TRUE)
})




# Applications of ICS for ICA ----

## Simulation of data----
X_ini <- simulation_ICA(n = Ntot)

## Computation of eigenvalues based on different methods ---
all_methods <- c("ICSEigen", "ICSQR")

lapply(alpha_all, function(alpha){

  eigenvalues_all <- sapply(k_all+1, function(k){
    eigenvalues_all_methods(X_ini, k, alpha)
  })
  eigenvalues_all_dt <- data.frame(eigenvalues_all)

  dt_eigenvalues_plot(all_methods, eigenvalues_all_dt, nb = ncol(X_ini),
                      alpha, df_name = "simulations_ICA_eigenvalues",
                      save = TRUE, k_all = k_all)
})



# Applications of ICS for outlier detection ----

## Nearly singular data set  ----
# Import data
df <- as.matrix(HTP3)
rownames(df) <- rownames(HTP3)
ind_outlier <- 32

print(format(kappa(df, exact = TRUE), scientific = TRUE, digits = 2))

# Boxplots
# description of scale of the absolute mean
df_range <- data.frame(range = apply(HTP3, 2 , function(x) abs(mean(x))))
df_range <- cbind(Variable = row.names(df_range), df_range)
df_range$interval = cut(df_range$range, breaks = 10^(seq(-15,15,3)),   
                        ordered_result = TRUE)
# table(df_range$interval)
HTP3_long <- reshape(HTP3, varying = list(1:ncol(HTP3)), v.names = "value",
               timevar = "Variable", direction = "long")
HTP3_long$Variable <- factor(paste0("V.", HTP3_long$Variable),
                             levels = colnames(HTP3), ordered = TRUE)
HTP3_long <- merge(HTP3_long, df_range, by = "Variable", sort = FALSE)

HTP3_box <- HTP3_long %>% ggplot(aes(x=Variable, y = value,
                                     color = interval)) +
  geom_boxplot() +
  coord_flip()+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1),
                     guide = guide_axis(angle = 45) ) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text.y = element_text(size = 10),
        panel.spacing = unit(1, "lines"))+
  scale_colour_colorblind()+
  facet_grid(~interval, scales = "free")


pdf("figures/HTP3_boxplot.pdf", width = 14)
print(HTP3_box)
dev.off()


# PCA
eigen(cov(scale(df, center = TRUE, scale = FALSE)), symmetric = TRUE, only.values = TRUE)


# ICS - ICS_EIGEN
res_ICS_EIGEN <- ICS::ics2(df)

# ICS - QR
sapply(alpha_all, function(alpha){
  res_ICS_QR <- ICS_QR_full_rank(Xt = df,
                                     alpha = alpha, cf = 1)
  # Eigenvalues

  plot.ics.eigenvalues(res_ICS_QR, df_name = "HTP3_nearly_singular_eigenvalues", save = TRUE)

  #eigen(cov(df))$values

  # Choice of number of components
  if(alpha < 0){
    k = ncol(res_ICS_QR$Z)
  }else{
    k = 1
  }

  plot.ics.distances(df, res_ICS_QR, ind_outlier, k,
                     df_name = "HTP3_nearly_singular_ICSD", save = TRUE)

})




## Collinear data set ----
# Import data
df <- as.matrix(HTP2)
rownames(df) <- rownames(HTP2)
ind_outlier <- 28

print(format(kappa(df, exact = TRUE), scientific = TRUE, digits = 2))

# Boxplots
# description of scale of the absolute mean
df_range <- data.frame(range = apply(HTP2, 2 , function(x) abs(mean(x))))
df_range <- cbind(Variable = row.names(df_range), df_range)
df_range$interval = cut(df_range$range, breaks = 10^(seq(-15,15,3)),              ordered_result = TRUE)
# table(df_range$interval)

HTP2_long <- reshape(HTP2, varying = list(1:ncol(HTP2)), v.names = "value",
                     timevar = "Variable", direction = "long")
HTP2_long$Variable <- factor(paste0("V.", HTP2_long$Variable),
                             levels = colnames(HTP2), ordered = TRUE)

HTP2_long <- merge(HTP2_long, df_range, by = "Variable", sort = FALSE)


HTP2_box <- HTP2_long %>% ggplot(aes(x=Variable, y = value,
                                     color = interval)) +
  geom_boxplot() +
  coord_flip()+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE,
                                                 digits = 1),
                     guide = guide_axis(angle = 45) ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text.y = element_text(size = 10),
        panel.spacing = unit(1, "lines"))+
  scale_colour_colorblind()+
  facet_grid(~interval, scales = "free")

pdf("figures/HTP2_boxplot.pdf", width = 14, height = 14)
print(HTP2_box)
dev.off()





# ICS - ICS_EIGEN

# PCA
eigen(cov(scale(df, center = TRUE, scale = FALSE)), symmetric = TRUE, only.values = TRUE)


# ICS - Eigen
# res_ICS_EIGEN <- ICS::ics2(df)
res_svd <- svd(df)
k <- sum((res_svd$d/res_svd$d[1]) > (max(dim(df)) *.Machine$double.eps))
X_reduc <- res_svd$u[,1:k]%*%diag(res_svd$d[1:k])
res_ICS_EIGEN <- ICS::ics2(X_reduc)

# ICS - QR
sapply(alpha_all, function(alpha){
  res_ICS_QR <- ICS_QR_not_full_rank(Xt = df, ref = "first",
                                     tol_eps = .Machine$double.eps,
                                     alpha = alpha, cf = 1)
  # Eigenvalues

  plot.ics.eigenvalues(res_ICS_QR, df_name = "HTP2_collinear_eigenvalues", save = TRUE)

  #eigen(cov(df))$values

  # Choice of number of components
  if(alpha < 0){
    k = ncol(res_ICS_QR$Z)-c(1:0)
  }else{
    k = 1
  }

  plot.ics.distances(df, res_ICS_QR, ind_outlier, k,
                     df_name = "HTP2_collinear_ICSD", save = TRUE)

})

