
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

inv_softmax <- function (x) {
  log(x) - max(log(x))
}


holmes_test = function(y, x, g, case_ind, ws=NULL){

  if(is.null(ws)){
    ws = rep(1, nrow(g))
  }

  y[!case_ind, ] = 0
  indicator = as.numeric(case_ind)
  data_all = data.frame(cbind(x,y,indicator,g))
  colnames(data_all) = c(colnames(x), colnames(y),'indicator','g')


  M0 = as.formula(paste('g~',
                        paste(colnames(x), collapse='+'),
                        '+1', sep=''))
  M1 = as.formula(paste('g~',
                        paste(colnames(x), collapse='+'),
                        '+indicator+1', sep=''))
  M2 = as.formula(paste('g~',
                        paste(colnames(x), collapse='+'),'+',
                        paste(colnames(y), collapse='+'),
                        '+indicator+1', sep=''))

  # fit models
  fit_M0 = glm(M0, data = data_all, family = 'binomial',weights = ws)
  fit_M1 = glm(M1, data = data_all, family = 'binomial',weights = ws)
  fit_M2 = glm(M2, data = data_all, family = 'binomial',weights = ws)

  Ftest01 = anova(fit_M0, fit_M1, test = "Chisq")
  Ftest12 = anova(fit_M1, fit_M2, test = "Chisq")
  Ftest02 = anova(fit_M0, fit_M2, test = "Chisq")

  out = list(pval01 = Ftest01$`Pr(>Chi)`[2],
             pval12 = Ftest12$`Pr(>Chi)`[2],
             pval02 = Ftest02$`Pr(>Chi)`[2])#,
             # coef_full_mod = summary(fit_M2)$coefficient)
  return(out)
}


holmes_test_simple = function(y, g, case_ind){

  y[!case_ind] = 0
  indicator = as.numeric(case_ind)
  data_all = data.frame(cbind(y,indicator,g))
  colnames(data_all) = c('y','indicator','g')


  M0 = as.formula('g~1')
  M1 = as.formula('g~indicator+1')
  M2 = as.formula('g~y+indicator+1')
  # fit models
  fit_M0 = glm(M0, data = data_all, family = 'binomial')
  fit_M1 = glm(M1, data = data_all, family = 'binomial')
  fit_M2 = glm(M2, data = data_all, family = 'binomial')

  Ftest01 = anova(fit_M0, fit_M1, test = "Chisq")
  Ftest02 = anova(fit_M0, fit_M2, test = "Chisq")
  Ftest12 = anova(fit_M1, fit_M2, test = "Chisq")

  out = c(Ftest01$`Pr(>Chi)`[2], Ftest02$`Pr(>Chi)`[2])
  # out = unlist(list(pval01 = Ftest01$`Pr(>Chi)`[2],
  #         pval02 = Ftest02$`Pr(>Chi)`[2],
  #         pval12 = Ftest12$`Pr(>Chi)`[2],
  #         M2coef = coef(fit_M2)['indicator'],
  #         M1coef = coef(fit_M1)['indicator']))
  if(is.na(out[1])) out[1] = 1
  return(out)
}

Np_fit = function(x, y, lambda = 1, BootStrap = T, BBs = 200){
  N = length(x)

  if(BootStrap){
    mwa = array(dim=c(BBs, N))
    for(bb in 1:BBs){
      ind = sample(1:N, N, T)
      xs = x[ind]
      ys = y[ind]
      xs = xs[order(ys)]
      ys = ys[order(ys)]
      for(i in 1:N){
        my_weights = exp(-(abs(ys[i] - ys)) * lambda)

        ws = my_weights/sum(my_weights)
        mwa[bb, i] = sum(xs*ws)
      }
    }
  } else {
    mwa = array(dim=N)
    for(i in 1:N){
      my_weights = exp(-(abs(y[i] - y)) * lambda)

      ws = my_weights/sum(my_weights)
      mwa[i] = sum(x*ws)
    }
  }
  return(mwa)
}

