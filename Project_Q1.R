setwd("S:/Acctn/Trinity/Modules/ASM/MainAssignment/Dataset")
library(plyr)
library(dplyr)
library(stats)
library(ggplot2)
library(MCMCpack)


###############################################
############################################### Q.1


data <- read.csv('winemag-data-130k-v2.csv')
View(data)

#### Below lines are used to plot points for few countries
country_data = data[(data$country == "US" | data$country == "Italy"| data$country == "France" | data$country == "Germany"
                     | data$country == 'South Africa' | data$country == "UK"),]

ggplot(country_data) + geom_boxplot(aes(country, points, fill = country))


### Loading data specific to Question 1
wine_data = data[((data$variety == "Sauvignon Blanc") & (data$country == 'South Africa') & (data$price ==15))
            |((data$country == 'Chile') & (data$variety == 'Chardonnay') &(data$price ==15)) ,]

View(wine_data)
count(wine_data, "country")
library(tidyr)
wine_data <- wine_data %>% drop_na()
count(wine_data, "variety")
wine_compare = wine_data[c("variety","points")]
wine_compare <- wine_compare %>% drop_na()

count(wine_compare,"variety")
ggplot(wine_compare) + geom_boxplot(aes(variety, points, fill = variety))+ geom_jitter(aes(variety, points, shape = wine_compare$variety))

mean.wine<-aggregate(points ~ variety, wine_compare, mean ) 
mean.sd<-aggregate(points ~ variety, wine_compare, sd ) 


t.test(points ~ variety, data=wine_compare, var.equal = TRUE)

mu =50
sd = 20
b0 = 100
a0 = 1

compare_2_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/(20^2), del0 = 0, gamma0 = 1/(20^2), a0 = 1, b0 = 100, maxiter = 5000)
{
  y1 <- y[ind == "Sauvignon Blanc"]
  y2 <- y[ind == "Chardonnay"]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

fit.data <- compare_2_gibbs(wine_compare$points, as.factor(wine_compare$variety))
plot(as.mcmc(fit.data ))
raftery.diag(as.mcmc(fit.data ))

apply(fit.data , 2, mean)
apply(fit.data , 2, sd)

# calculating probablilty 
y11_sim <- rnorm(5000, fit.data[, 1] + fit.data[, 2], sd = 1/sqrt(fit.data[, 3]))
y12_sim <- rnorm(5000, fit.data[, 1] - fit.data[, 2], sd = 1/sqrt(fit.data[, 3]))
ggplot(data.frame(sim_diff = y11_sim - y12_sim), aes(x=sim_diff)) + stat_bin(aes(sim_diff)) +geom_histogram(color="blue", fill="white")

mean(y11_sim > y12_sim)
ggplot(data.frame(y11_sim, y12_sim)) + geom_point(color='#8dd3c7',fill="white",
                                                  aes(y11_sim, y12_sim), alpha = 0.3) + geom_abline(slope = 1, intercept = 0) 

# how much better?
(mean(y11_sim) - mean(y12_sim))

## Probability that Sauvignon Blanc will be better
mean(y11_sim > y12_sim)


##################################################
###########################################Question 2
italy_wine <- data[(data$country == 'Italy') & (data$price < 20),]
View(italy_wine)

italy_wine=italy_wine[rowSums(is.na(italy_wine)) == 0,]

italy_wine <- italy_wine[c('points','title','region_1','variety')]
##Picking those regions which has at least 4 reviews
x <-italy_wine %>% group_by(region_1) %>% summarise(four_reviews = n_distinct(title)) %>% filter(four_reviews >=4)
summary(x)
drop_x<-x$region_1==""
x<-x[!drop_x,]
x<-as.data.frame(x)

## Selecting those regions with wines costing than 20 and having at least 4 reviews
selected_region <-italy_wine[italy_wine$region_1 %in% x$region_1,]

selected_region$ind = as.numeric(factor(selected_region$region_1))

selected_region<-selected_region[order(selected_region$ind),]


ggplot(selected_region) + geom_boxplot(aes(x = reorder(region_1, points, median), points,
                               fill = reorder(region_1, points, median)), show.legend=FALSE)

ggplot(selected_region, aes(x = reorder(region_1, region_1, length))) + stat_count()
ggplot(selected_region, aes(points)) + stat_bin()

compare_m_gibbs2 <- function(y, ind, mu0 = 50, tau0 = 1/400,
                            a0 = 1, b0 = 50, alpha0 =1, beta0 = 50, maxiter = 5000)
{
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 100 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  ###
  ### starting values
  m <- nlevels(ind)
  ### Changed this line to 2 as only ywo countries are selected
  
  ybar <- theta <- tapply(y, ind, mean)
  x <- tapply(y, ind, var)
  x <- x[!is.na(x)]
  ## This block is there to handle na and infinite records
  tau_w <- mean(1/x[!is.infinite(1/x)]) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  alphan <- alpha0 + sum(n_m)/2
  ###
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  ### MCMC algorithm
  for(s in 1:maxiter)
  {
    # sample new values of the thetas
    for(j in 1:m)
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    betan <- beta0 + ss/2
    tau_w <- rgamma(1, alphan, betan)
    #sample a new value of mu
    taum <- m * tau_b + tau0
    mum <- (mean(theta) * m * tau_b + mu0 * tau0) / taum
    mu <- rnorm(1, mum, 1/ sqrt(taum))
    # sample a new value of tau_b
    am <- a0 + m/2
    bm <- b0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, am, bm)
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

selected_region$ind<-as.factor(selected_region$ind)
fit2.data <- compare_m_gibbs2(selected_region$points, selected_region$ind)
apply(fit2.data$params, 2, mean)
apply(fit2.data$params, 2, sd)

## within regions standard variation
mean(1/sqrt(fit2.data$params[, 2]))
sd(1/sqrt(fit2.data$params[, 2]))

## between regions standard variation
mean(1/sqrt(fit2.data$params[, 3]))
sd(1/sqrt(fit2.data$params[, 3]))

theta_hat <- apply(fit2.data$theta, 2, mean) ## posterior mean summary among 152 regions 
#names(theta_hat)<-1:152
sort(theta_hat, decreasing = TRUE)

## Now to check regions which has ratings better than average wine
Above_Avg<-which(theta_hat>mean(theta_hat))
Above_AvgRegion <-unique(selected_region[selected_region$ind %in% Above_Avg,]$region_1)
print(Above_AvgRegion)
## Total 73 regions above average rating

## Check for region which has the highest rating.
which(theta_hat == max(theta_hat))
selected_region[selected_region$ind == 126,]
## Trento is the region having maximum rating value


theta_ci <- apply(fit2.data$theta, 2, quantile, prob = c(0.025, .975))
df_error <- data.frame(lower = theta_ci[1, ], upper = theta_ci[2, ], mean = theta_hat,
                       regions = factor(1:152))



ggplot(df_error, aes(x = reorder(regions, mean), mean)) + geom_errorbar(aes(ymin = lower, ymax = upper))

theta_df <- data.frame(samples = as.numeric(fit2.data$theta),
                       regions = rep(1:ncol(fit2.data$theta), each = nrow(fit2.data$theta)))
ggplot(theta_df) + geom_boxplot(aes(x = reorder(regions, samples, median), samples,
                                    fill = reorder(regions, samples, median)), show.legend=FALSE)


y_bar= tapply(selected_region$points, selected_region$ind, mean)-theta_hat
ggplot(data.frame(size = tapply(selected_region$points, selected_region$ind, length), y_bar = y_bar), aes(size, y_bar)) + geom_point()



ggplot(data.frame(size = tapply(selected_region$points, selected_region$ind, length), theta_hat = theta_hat),
       aes(size, theta_hat)) + geom_point()


ggplot(data.frame(ybar = tapply(selected_region$points, selected_region$ind, mean), theta_hat = theta_hat),
       aes(ybar, theta_hat)) + geom_point()

################################################



