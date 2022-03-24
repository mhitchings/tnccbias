### Bias correction for test-negative case-control studies

# Function to calculate VE without any bias correction
VE_biased = function(t,phi,theta,alpha_vi,alpha_ui,alpha_vn,alpha_un,lambda_i) {
  1 - alpha_un/alpha_vn * 
    ((1-phi)*(1-exp(-alpha_vi*lambda_i*t)) + phi*(1-exp(-alpha_vi*theta*lambda_i*t)))/((1-exp(-alpha_ui*lambda_i*t)))
}

# Function to calculate bias-corrected VE (Equation 2)
VE_biascorrected = function(t,Tp,Tv,phi,theta,alpha_vfi,alpha_ui,alpha_vfn,alpha_un,lambda_i,
                               phi_tinit,theta_tinit,
                               alpha_vri,alpha_vrn) {
  
  1-(1-VE_biased(t-Tp,phi,theta,alpha_vfi,alpha_ui,alpha_vfn,alpha_un,lambda_i))/
    (1-VE_biased(Tp-Tv,phi_tinit,theta_tinit,alpha_vri,alpha_ui,alpha_vrn,alpha_un,lambda_i))
}

# Function to calculate the bias indicator
biasindicator = function(t,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                         phi_tinit,theta_tinit) {
  cri = v*((1-phi_tinit)*(1-exp(-alpha_ri*lambda_i*Tp)) + phi_tinit*(1-exp(-alpha_ri*lambda_i*theta_tinit*Tp)))
  crn = v*alpha_rn*Tp
  cui = (1-v)*(1-exp(-alpha_ui*lambda_i*t)) + v*(1-exp(-alpha_pi*lambda_i*Tv))
  cun = (1-v)*alpha_un*t + v*alpha_pn*Tv
  
  return(cun*cri/cui/crn)
}

# Function to calculate bias indicator when there is differential symptom probability or testing probability
# between vaccinated groups
biasindicator_testing = function(t,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                 phi_tinit,theta_tinit,
                                 pi_rim,pi_ris,pi_rnm,pi_rns,pi_uim,pi_uis,pi_unm,pi_uns,
                                 mu_rm,mu_rs,mu_um,mu_us,
                                 xi_m,xi_s) {
  
  return(biasindicator(t,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                       phi_tinit,theta_tinit) * 
           (pi_rim*mu_rm*xi_m+pi_ris*mu_rs*xi_s) * (pi_unm*mu_um*xi_m+pi_uns*mu_us*xi_s) / (pi_rnm*mu_rm*xi_m+pi_rns*mu_rs*xi_s) / (pi_uim*mu_um*xi_m+pi_uis*mu_us*xi_s))
  
}

# Function to calculate VE without bias correction when there is differential symptom probability or testing probability
# between vaccinated groups
VE_biased_testing = function(t,phi,theta,alpha_vi,alpha_ui,alpha_vn,alpha_un,lambda_i,
                     pi_vim,pi_vis,pi_vnm,pi_vns,pi_uim,pi_uis,pi_unm,pi_uns,
                     mu_vm,mu_vs,mu_um,mu_us,
                     xi_m,xi_s) {
  1 - alpha_un/alpha_vn * 
    ((1-phi)*(1-exp(-alpha_vi*lambda_i*t)) + phi*(1-exp(-alpha_vi*theta*lambda_i*t)))/((1-exp(-alpha_ui*lambda_i*t))) * 
    (pi_vim*mu_vm*xi_m+pi_vis*mu_vs*xi_s) * (pi_unm*mu_um*xi_m+pi_uns*mu_us*xi_s) / (pi_vnm*mu_vm*xi_m+pi_vns*mu_vs*xi_s) / (pi_uim*mu_um*xi_m+pi_uis*mu_us*xi_s)
}

# Function to calculate bias-corrected VE when there is differential symptom probability or testing probability
# between vaccinated groups
VE_biascorrected_testing = function(t,Tp,Tv,phi,theta,alpha_vfi,alpha_ui,alpha_vfn,alpha_un,lambda_i,
                               phi_tinit,theta_tinit,
                               alpha_vri,alpha_vrn,
                               pi_vfim,pi_vfis,pi_vfnm,pi_vfns,
                               pi_vrim,pi_vris,pi_vrnm,pi_vrns,
                               pi_uim,pi_uis,pi_unm,pi_uns,
                               mu_vfm,mu_vfs,
                               mu_vrm,mu_vrs,
                               mu_um,mu_us,
                               xi_m,xi_s) {
  
  1-(1-VE_biased_testing(t-Tp,phi,theta,alpha_vfi,alpha_ui,alpha_vfn,alpha_un,lambda_i,
                 pi_vfim,pi_vfis,pi_vfnm,pi_vfns,pi_uim,pi_uis,pi_unm,pi_uns,
                 mu_vfm,mu_vfs,mu_um,mu_us,
                 xi_m,xi_s))/
    (1-VE_biased_testing(Tp-Tv,phi_tinit,theta_tinit,alpha_vri,alpha_ui,alpha_vrn,alpha_un,lambda_i,
                 pi_vrim,pi_vris,pi_vrnm,pi_vrns,pi_uim,pi_uis,pi_unm,pi_uns,
                 mu_vrm,mu_vrs,mu_um,mu_us,
                 xi_m,xi_s))
}


# Plots
par(mar=c(4,4,1,1))

ylim_bi=c(0.9,1.7)
ylim_bc=c(0.5,1)
# Left column: bias indicator
# Right column: bias-correction

# Row 1: time-invariant bias
# Row 2: pending vaccination have lower risk
# Row 3: pending vaccination and recently vaccinated have lower risk
# Row 4: vaccine has some effect in initial period
# Row 5: early vs. late adopters?

# Parameters
phi=1
theta=0.3
v=0.7

phi_tinit=0
theta_tinit=1

alpha_fi = 1.25
alpha_ui = 1
alpha_fn = 1
alpha_un = 1

alpha_ri=1.25
alpha_rn=1

alpha_pi=1.25
alpha_pn=1

lambda_i = 0.005
Tv=Tve=30
Tp=7

ts_toplot = (Tv+Tp+1):100

par(mfrow=c(5,2))
# Time-invariant alpha: higher risk among fully, recently, and pending vaccinated
plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) biasindicator(x+Tp+Tv,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                                                           phi_tinit,theta_tinit)),
     type='l',ylim=ylim_bi,ylab='Bias indicator',xlab='t')
abline(h=1,lwd=2,lty=2)


plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biased(x+Tp+Tv,phi,theta,alpha_fi,alpha_ui,alpha_fn,alpha_un,lambda_i)),
     type='l',ylim=ylim_bc,ylab='VE',xlab='t',lty=1)
abline(h=1-theta,lwd=2,lty=2)
lines(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biascorrected(x+Tp+Tv,Tve+Tp,Tve,phi,theta,alpha_ri,alpha_ui,alpha_rn,alpha_un,lambda_i,
                                                                 phi_tinit,theta_tinit,
                                                                 alpha_fi,alpha_fn)),
      lty=3)

# Time-varying alpha: higher risk in recently vaccinated only
alpha_pi=alpha_fi=1

plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) biasindicator(x+Tp+Tv,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                                                      phi_tinit,theta_tinit)),
     type='l',ylim=ylim_bi,ylab='Bias indicator',xlab='t')
abline(h=1,lwd=2,lty=2)


plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biased(x+Tp+Tv,phi,theta,alpha_fi,alpha_ui,alpha_fn,alpha_un,lambda_i)),
     type='l',ylim=ylim_bc,ylab='VE',xlab='t',lty=1)
abline(h=1-theta,lwd=2,lty=2)
lines(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biascorrected(x+Tp+Tv,Tve+Tp,Tve,phi,theta,alpha_ri,alpha_ui,alpha_rn,alpha_un,lambda_i,
                                                                             phi_tinit,theta_tinit,
                                                                             alpha_fi,alpha_fn)),
      lty=3)

# Time-varying alpha: lower risk among recently vaccinated and pending vaccination than fully vaccinated
alpha_ri=0.8
alpha_pi=0.8

plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) biasindicator(x+Tp+Tv,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                                          phi_tinit,theta_tinit)),
     type='l',ylim=ylim_bi,ylab='Bias indicator',xlab='t')
abline(h=1,lwd=2,lty=2)


plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biased(x+Tp+Tv,phi,theta,alpha_fi,alpha_ui,alpha_fn,alpha_un,lambda_i)),
     type='l',ylim=ylim_bc,ylab='VE',xlab='t',lty=1)
abline(h=1-theta,lwd=2,lty=2)
lines(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biascorrected(x+Tp+Tv,Tve+Tp,Tve,phi,theta,alpha_ri,alpha_ui,alpha_rn,alpha_un,lambda_i,
                                                                 phi_tinit,theta_tinit,
                                                                 alpha_fi,alpha_fn)),
      lty=3)

# Some vaccine effectiveness in initial period
alpha_ri=1.25
alpha_pi=1.25
alpha_fi=1.25
theta_tinit=0.9
phi_tinit=1

plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) biasindicator(x+Tp+Tv,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                                          phi_tinit,theta_tinit)),
     type='l',ylim=ylim_bi,ylab='Bias indicator',xlab='t')
abline(h=1,lwd=2,lty=2)


plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biased(x+Tp+Tv,phi,theta,alpha_fi,alpha_ui,alpha_fn,alpha_un,lambda_i)),
     type='l',ylim=ylim_bc,ylab='VE',xlab='t',lty=1)
abline(h=1-theta,lwd=2,lty=2)
lines(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biascorrected(x+Tp+Tv,Tve+Tp,Tve,phi,theta,alpha_ri,alpha_ui,alpha_rn,alpha_un,lambda_i,
                                                                 phi_tinit,theta_tinit,
                                                                 alpha_fi,alpha_fn)),
      lty=3)

# Time-varying changes in testing behaviour: no testing for moderate symptoms in recently/fully vaccinated
theta_tinit=1
phi_tinit=0

pi_fim=0.5
pi_fis=0.04
pi_fnm=0.5
pi_fns=0.025
pi_rim=0.5
pi_ris=0.04
pi_rnm=0.5
pi_rns=0.025
pi_uim=0.5
pi_uis=0.04
pi_unm=0.5
pi_uns=0.025

mu_fm=0.025
mu_fs=0.1
mu_rm=0.025
mu_rs=0.1
mu_um=0.1
mu_us=0.1

xi_m=0.5
xi_s=1

plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) biasindicator_testing(x+Tp+Tv,Tp,Tv,phi,theta,alpha_ri,alpha_pi,alpha_ui,alpha_rn,alpha_pn,alpha_un,lambda_i,v,
                                                                  phi_tinit,theta_tinit,
                                                                  pi_rim,pi_ris,pi_rnm,pi_rns,pi_uim,pi_uis,pi_unm,pi_uns,
                                                                  mu_rm,mu_rs,mu_um,mu_us,
                                                                  xi_m,xi_s)),
     type='l',ylim=ylim_bi,ylab='Bias indicator',xlab='t')
abline(h=1,lwd=2,lty=2)

plot(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biased_testing(x+Tp+Tv,phi,theta,alpha_fi,alpha_ui,alpha_fn,alpha_un,lambda_i,
                                                                    pi_fim,pi_fis,pi_fnm,pi_fns,pi_uim,pi_uis,pi_unm,pi_uns,
                                                                    mu_fm,mu_fs,mu_um,mu_us,
                                                                    xi_m,xi_s)),
     type='l',ylim=ylim_bc,ylab='VE',xlab='t',lty=1)
abline(h=1-theta,lwd=2,lty=2)
lines(ts_toplot-Tv-Tp,sapply(ts_toplot-Tv-Tp,function(x) VE_biascorrected_testing(x+Tp+Tv,Tve+Tp,Tve,phi,theta,alpha_ri,alpha_ui,alpha_rn,alpha_un,lambda_i,
                                                                 phi_tinit,theta_tinit,
                                                                 alpha_fi,alpha_fn,
                                                                 pi_fim,pi_fis,pi_fnm,pi_fns,
                                                                 pi_rim,pi_ris,pi_rnm,pi_rns,
                                                                 pi_uim,pi_uis,pi_unm,pi_uns,
                                                                 mu_fm,mu_fs,
                                                                 mu_rm,mu_rs,
                                                                 mu_um,mu_us,
                                                                 xi_m,xi_s)),
      lty=3)

