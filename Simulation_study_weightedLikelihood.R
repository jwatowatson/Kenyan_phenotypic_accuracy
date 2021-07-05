### Simulations to show how re-weighted likelihoods work ####
rm(list=ls())
allele_freq = 0.1 # allele frequency in controls
p_effect = 0.7 # reduction of allele frequency in true cases

true_logodds = log(allele_freq/(1-allele_freq)) - log((p_effect*allele_freq)/(1-p_effect*allele_freq))

prop_contam = seq(0.5, 1, length.out = 20)
my_beta = 0.2
N = 1000
Nsims = 1000

res_weighted = res_standard = array(dim = c(length(prop_contam),2, Nsims))
pb=txtProgressBar(min=0,max = length(prop_contam), style = 3)
type2error = array(dim = c(length(prop_contam),2, Nsims))
for(i in 1:length(prop_contam)){
  N1 = round(prop_contam[i]*N)
  N2 = N-N1
  for(k in 1:Nsims){
    g_controls = rbinom(N, size=2, prob = allele_freq)
    my_ws = c(rbeta(n = N1, shape1 = 1, shape2 = my_beta),
              rbeta(n = N2, shape1 = my_beta, shape2 = 1))

    true_state = c(rep(1, N1), rep(0, N2))
    allele_freq_cont = ifelse(true_state==1, allele_freq*p_effect, allele_freq)
    g_cases = rbinom(N, size=2, prob = allele_freq_cont)

    mod_w = glm(cc ~ g, family = 'quasibinomial',
                data = data.frame(cc=c(rep(0,N),rep(1,N)),
                                  g = c(g_controls, g_cases)),
                weights = c(rep(1, N), my_ws))

    # compute the standard error for the weighted model - takes into account reduced effective sample size
    cov.m1 <- sandwich::vcovHC(mod_w, type = "HC0")
    std.err <- sqrt(diag(cov.m1))['g']
    res_weighted[i,1:2,k] = c(coef(mod_w)['g'], std.err)
    type2error[i,1,k] = as.numeric(sign(coef(mod_w)['g'] + 1.96*std.err) == sign(coef(mod_w)['g'] - 1.96*std.err))


    mod_stand = glm(cc ~ g, family = 'quasibinomial',
                    data = data.frame(cc=c(rep(0,N),rep(1,N)),
                                      g = c(g_controls, g_cases)))
    res_standard[i,1:2,k] = summary(mod_stand)$coefficients['g',1:2]
    std.err = summary(mod_stand)$coefficients['g',2]
    type2error[i,2,k] = as.numeric(sign(coef(mod_stand)['g'] + 1.96*std.err) == sign(coef(mod_stand)['g'] - 1.96*std.err))

  }
  setTxtProgressBar(pb, i)
}


xs=100*prop_contam 
pdf('sim_study_results1.pdf', height = 10, width = 10)
par(mfrow=c(2,2), las=1, cex.lab=1.3, cex.axis=1.1,mar=c(5,7,4,2))

hist(c(rbeta(n = 500, shape1 = 1, shape2 = my_beta),
       rbeta(n = 500, shape1 = my_beta, shape2 = 1)),
     main='', xlab = 'Case probability weights',breaks = 20)
mtext(side = 3, text = 'A',line = 2, adj = 0, cex = 1.5)
hist(c(rbeta(n = 1000, shape1 = 1, shape2 = my_beta),
       rbeta(n = 0, shape1 = my_beta, shape2 = 1)),
     main='', xlab = 'Case probability weights', breaks = 20)
mtext(side = 3, text = 'B',line = 2, adj = 0, cex = 1.5)

plot(xs, abs(rowMeans(res_standard[,1,])),
     ylim=c(0.15,0.4),panel.first = grid(),
     xlab='Proportion of true cases (%)',
     ylab = '', pch=16,xlim=c(50,100))
mtext(text = 'Mean absolute effect size (log-odds)',
      side = 2,line = 5,las=3, cex = 1.1)
mtext(side = 3, text = 'C',line = 2, adj = 0, cex = 1.5)

abline(h = true_logodds)
points(xs, abs(rowMeans(res_weighted[,1,])),
       pch=17,col='red')

plot(xs, rowMeans(res_standard[,2,]),
     ylim = c(0.109, .135),pch=16,
     xlab='Proportion of true cases (%)',
     ylab = '',panel.first = grid())
mtext(text = 'Mean standard error (log-odds scale)',
      side = 2,line = 5,las=3, cex = 1.1)
points(xs, rowMeans(res_weighted[,2,]),col='red',pch=17)
mtext(side = 3, text = 'D',line = 2, adj = 0, cex = 1.5)

legend('topright',pch=16:17, legend = c('Non-weighted','Weighted'),
       col = c(1,2), inset=0.03)
dev.off()



pdf('sim_study_results2.pdf', height = 8, width = 8)
par(mfrow=c(1,1), las=1, cex.lab=1.3, cex.axis=1.3,mar=c(5,7,2,2))
plot(xs, 100*rowMeans(type2error[,1,]), ylim = c(0,100),type='l',lwd=3,
     ylab = '1 - Type 2 error (power, %)',panel.first = grid(),
     xlab='Proportion of true cases (%)',col='red')
lines(xs, 100*rowMeans(type2error[,2,]), lty=2, lwd=3)
dev.off()
