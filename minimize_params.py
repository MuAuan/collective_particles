#https://qiita.com/Student-M/items/4e3e286bf08b7320b665

#include package
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#define differencial equation of seir model
def seir_eq(v,t,beta,lp,ip):
    N=763
    a = -beta*v[0]*v[2]/N
    b = beta*v[0]*v[2]/N-(1/lp)*v[1]
    c = (1/lp)*v[1]-(1/ip)*v[2]
    d = (1/ip)*v[2]
    return [a,b,c,d]

#solve seir model
N,S0,E0,I0=762,0,1,0
ini_state=[N,S0,E0,I0]
beta,lp,ip=6.87636378, 1.21965986, 2.01373496 #2.493913  , 0.95107715, 1.55007883
t_max=14
dt=0.01
t=np.arange(0,t_max,dt)
plt.plot(t,odeint(seir_eq,ini_state,t,args=(beta,lp,ip))) #0.0001,1,3
plt.legend(['Susceptible','Exposed','Infected','Recovered'])
plt.pause(1)
plt.close()

#show observed i
#obs_i=np.loadtxt('fitting.csv')
data_influ=[3,8,28,75,221,291,255,235,190,125,70,28,12,5]
data_day = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
obs_i = data_influ

plt.plot(obs_i,"o", color="red",label = "data")
plt.legend()
plt.pause(1)
plt.close()

#function which estimate i from seir model func 
def estimate_i(ini_state,beta,lp,ip):
    v=odeint(seir_eq,ini_state,t,args=(beta,lp,ip))
    est=v[0:int(t_max/dt):int(1/dt)]
    return est[:,2]
    
#define logscale likelihood function
def y(params):
    est_i=estimate_i(ini_state,params[0],params[1],params[2])
    return np.sum(est_i-obs_i*np.log(np.abs(est_i)))
    
#optimize logscale likelihood function
mnmz=minimize(y,[beta,lp,ip],method="nelder-mead")
print(mnmz)
#R0
#N_total = S_0+I_0+R_0
#R0 = N_total*beta_const *(1/gamma_const)
beta_const,lp,gamma_const = mnmz.x[0],mnmz.x[1],mnmz.x[2] #感染率、感染待時間、除去率（回復率）
print(beta_const,lp,gamma_const)
R0 = beta_const*(1/gamma_const)
print(R0)

#plot reult with observed data
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
lns1=ax1.plot(obs_i,"o", color="red",label = "data")
lns2=ax1.plot(estimate_i(ini_state,mnmz.x[0],mnmz.x[1],mnmz.x[2]), label = "estimation")
lns_ax1 = lns1+lns2
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)

lns3=ax2.plot(obs_i,"o", color="red",label = "data")
lns4=ax2.plot(t,odeint(seir_eq,ini_state,t,args=(mnmz.x[0],mnmz.x[1],mnmz.x[2])))
ax2.legend(['data','Susceptible','Exposed','Infected','Recovered'], loc=0)
ax2.set_title('SEIR_b{:.2f}_ip{:.2f}_gamma{:.2f}_N{:d}_E0{:d}_I0{:d}_R0{:.2f}'.format(beta_const,lp,gamma_const,N,E0,I0,R0))
plt.savefig('./fig/SEIR_b{:.2f}_ip{:.2f}_gamma{:.2f}_N{:d}_E0{:d}_I0{:d}_R0{:.2f}_.png'.format(beta_const,lp,gamma_const,N,E0,I0,R0)) 
plt.show()
plt.close()