first_global_DOC<-function(data){

  #creating function for parameters estimation
  estim<-function(fr){
    #defining the objective function
    est<-function(x){

      parnames<-c("ks", "Eas")
      par<-x[1:length(parnames)]

      Ctot<-as.numeric(fr[, "DOC"])
      Temp<-as.numeric(fr[, "Temperature"])

      yhat<-Ctot*par[1]*exp(par[2]/Temp/8.314)
      obs<-as.numeric(fr[, "Rsoil"])

      # SSres<-sum(((obs-yhat)^2), na.rm = T)
      # SStot<-sum(((obs-mean(obs, na.rm=T))^2), na.rm = T)
      # Rsq<-1-(SSres/SStot)
      ll<--length(obs)*log(2*mean(obs, na.rm=T)*sd(obs, na.rm=T)^2)/2-sum((obs-yhat)^2, na.rm=T)/2/sd(obs, na.rm=T)^2
      # AIC<-2*length(x)-2*ll

      return(-2*ll)

    }

    #defining the goodness of fit function
    goodness<-function(x){

      parnames<-c("k", "Eas")
      par<-x[1:length(parnames)]

      Ctot<-as.numeric(fr[, "DOC"])
      Temp<-as.numeric(fr[, "Temperature"])

      yhat<-Ctot*par[1]*exp(par[2]/Temp/8.314)
      obs<-as.numeric(fr[, "Rsoil"])

      SSres<-sum(((obs-yhat)^2), na.rm = T)
      SStot<-sum(((obs-mean(obs, na.rm=T))^2), na.rm = T)
      Rsq<-1-(SSres/SStot)
      ll<--length(obs)*log(2*mean(obs, na.rm=T)*sd(obs, na.rm=T)^2)/2-sum((obs-yhat)^2, na.rm=T)/2/sd(obs, na.rm=T)^2
      AIC<-2*length(x)-2*ll

      return(c(R2=Rsq,AIC=AIC, ll=ll))

    }

    #approximate parameter estimation is done by MCMC method
    par_mcmc<-modMCMC(f=est, p=c(k=1e10, Eas=-50000),
                      lower=c(k=1e-6, Eas=-3000000),
                      upper=c(k=1e50, Eas=10000), niter=50000)

    #lower and upper limits for parameters are extracted
    pl<-summary(par_mcmc)["min",]
    pu<-summary(par_mcmc)["max",]

    #these limits are used to find global optimum by DEoptim
    opt_par<-DEoptim(fn=est, lower=pl, upper=pu,
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8,
                                 trace=FALSE, strategy=3, NP=500))

    #goodness of fit
    fit<-goodness(opt_par$optim$bestmem)

    #approximate parameter estimation is done by MCMC method
    par_prof<-modMCMC(f=est, p=opt_par$optim$bestmem,
                      lower=pl,
                      upper=pu, niter=50000)

    #goodness of fit
    fit<-goodness(opt_par$optim$bestmem)

    #best parameters
    p<-opt_par$optim$bestmem
    names(p)<-c("ks", "Eas")

    #sd of parameters
    p.sd<-summary(par_prof)[2,]

    #return list with opt_par and par_prof
    estim_out<-list()
    estim_out$pars<-p
    estim_out$pars.sd<-p.sd
    estim_out$fit<-fit

    return(estim_out)
  }

  #do the calculation
  res<-estim(fr=data)
  ks<-as.numeric(res$pars[1])
  Eas<-as.numeric(res$pars[2])
  ks.sd<-as.numeric(res$pars.sd[1])
  Eas.sd<-as.numeric(res$pars.sd[2])
  R2s<-as.numeric(res$fit[1])
  AICs<-as.numeric(res$fit[2])
  lls<-as.numeric(res$fit[3])



  return(c(ks=ks, Eas=Eas, ks.sd=ks.sd, Eas.sd=Eas.sd, R2=R2s, AIC=AICs, ll=lls))

}
