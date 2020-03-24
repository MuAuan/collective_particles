#https://qiita.com/Student-M/items/4e3e286bf08b7320b665

#include package
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import pandas as pd

#pandasでCSVデータ読む。
data = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv')
data_r = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv')
data_d = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv')

    
confirmed = [0] * (len(data.columns) - 4)
confirmed_r = [0] * (len(data_r.columns) - 4)
confirmed_d = [0] * (len(data_d.columns) - 4)
days_from_22_Jan_20 = np.arange(0, len(data.columns) - 4, 1)

#City
#city = "Hubei"
#city = "Korea, South"
#city = "Italy"
#city = "Spain"
city = "Iran"
#city = "Japan"

#データを加工する
t_cases = 0

for i in range(0, len(data), 1):
    if (data.iloc[i][1] == city): #for country/region
    #if (data.iloc[i][0] == city):  #for province:/state  
        print(str(data.iloc[i][0]) + " of " + data.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed[day - 4] += data.iloc[i][day] -  data_r.iloc[i][day] #+=
            confirmed_r[day - 4] += data_r.iloc[i][day]
            confirmed_d[day - 4] += data_d.iloc[i][day]*1 #fro drawings

#define differencial equation of seir model
def seir_eq(v,t,beta,lp,ip,N0):
    N=N0 #int(26749*1) #58403*4
    a = -beta*v[0]*v[2]/N
    b = beta*v[0]*v[2]/N-(1/lp)*v[1] #N
    c = (1/lp)*v[1]-(1/ip)*v[2]
    d = (1/ip)*v[2]
    return [a,b,c,d]

#solve seir model
N0 = 22116 #1042 japan
N,S0,E0,I0=int(N0*1.8),0,int(1),int(0)
ini_state=[N,S0,E0,I0]
beta,lp,ip, N= 0.10062997064586765, 1.0996353042201308e-13, 20.613246427653134, 14330.117810270265 #0.07517828536584859, 5.902996787735029e-14, 271.34861852988956 #japan #0.0349575134313698, 8.229631546593396e-10, 157.172343954462 #0.2,0.001,200 #1,2,7.4
t_max=60 #len(days_from_22_Jan_20)
dt=0.01
t=np.arange(0,t_max,dt)
plt.plot(t,odeint(seir_eq,ini_state,t,args=(beta,lp,ip,N))) #0.0001,1,3
plt.legend(['Susceptible','Exposed','Infected','Recovered'])
plt.pause(1)
plt.close()

#show observed i
#obs_i=np.loadtxt('fitting.csv')
#data_influ=[3,8,28,75,221,291,255,235,190,125,70,28,12,5]
#data_day = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#obs_i = data_influ
obs_i = confirmed

plt.plot(obs_i,"o", color="red",label = "data")
plt.legend()
plt.pause(1)
plt.close()

#function which estimate i from seir model func 
def estimate_i(ini_state,beta,lp,ip,N):
    v=odeint(seir_eq,ini_state,t,args=(beta,lp,ip,N))
    est=v[0:int(t_max/dt):int(1/dt)]
    return est[:,2]
    
#define logscale likelihood function
def y(params):
    est_i=estimate_i(ini_state,params[0],params[1],params[2],params[3])
    return np.sum(est_i-obs_i*np.log(np.abs(est_i)))
    
#optimize logscale likelihood function
mnmz=minimize(y,[beta,lp,ip,N],method="nelder-mead")
print(mnmz)
#R0
#N_total = S_0+I_0+R_0
#R0 = N_total*beta_const *(1/gamma_const)
beta_const,lp,gamma_const,N = mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3] #感染率、感染待時間、除去率（回復率）
print(beta_const,lp,gamma_const,N)
R0 = N*beta_const*(1/gamma_const)
print(R0)

#plot reult with observed data
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
lns1=ax1.plot(obs_i,"o", color="red",label = "data")
lns2=ax1.plot(confirmed_r,"*", color="green",label = "data_r")
lns3=ax1.plot(estimate_i(ini_state,mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3]), label = "estimation")
lns0=ax1.plot(t,odeint(seir_eq,ini_state,t,args=(mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3])))
lns_ax1 = lns1+lns2 + lns3 + lns0
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)
ax1.set_ylim(0,N0)

t_max=200 #len(days_from_22_Jan_20)
t=np.arange(0,t_max,dt)

lns4=ax2.plot(obs_i,"o", color="red",label = "data")
lns5=ax2.plot(confirmed_r,"*", color="green",label = "recovered")
lns6=ax2.plot(t,odeint(seir_eq,ini_state,t,args=(mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3])))
ax2.legend(['data','data_r','Susceptible','Exposed','Infected','Recovered'], loc=1)
ax2.set_title('SEIR_b{:.2f}_ip{:.2f}_gamma{:.2f}_N{:.2f}_E0{:d}_I0{:d}_R0{:.2f}'.format(beta_const,lp,gamma_const,N,E0,I0,R0))
ax1.grid()
ax2.grid()
plt.savefig('./fig/SEIR_{}_b{:.2f}_ip{:.2f}_gamma{:.2f}_N{:.2f}_E0{:d}_I0{:d}_R0{:.2f}_.png'.format(city,beta_const,lp,gamma_const,N,E0,I0,R0)) 
plt.show()
plt.close()