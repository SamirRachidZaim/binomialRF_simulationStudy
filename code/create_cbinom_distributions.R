#### Create cbinom distributions


### 500 trials 
cbinom <- correlbinom::correlbinom(.25, successprob = 1/20000, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_20000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/10000, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_10000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/5000, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_5000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/2000, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_2000F_.25Rho.csv')


cbinom <- correlbinom::correlbinom(.25, successprob = 1/1000, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_1000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/500, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_500F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/100, trials = 500, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/500T_100F_.25Rho.csv')







### 1000 trials 
cbinom <- correlbinom::correlbinom(.25, successprob = 1/20000, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_20000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/10000, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_10000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/5000, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_5000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/2000, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_2000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/1000, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_1000F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/500, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_500F_.25Rho.csv')

cbinom <- correlbinom::correlbinom(.25, successprob = 1/100, trials = 1000, model = 'witt')
write.csv(cbinom, 'cbinom_distributions/1000T_100F_.25Rho.csv')
