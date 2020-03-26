#https://qiita.com/Student-M/items/4e3e286bf08b7320b665

#include package
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import pandas as pd

#pandasでCSVデータ読む。
data = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv')
data_r = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv')
data_d = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv')

    
confirmed = [0] * (len(data.columns) - 4)
confirmed_r = [0] * (len(data_r.columns) - 4)
confirmed_d = [0] * (len(data_d.columns) - 4)
days_from_22_Jan_20 = np.arange(0, len(data.columns) - 4, 1)

#City
city = "Hubei"
#city = "Korea, South"
#city = "Italy"
#city = "Spain"
#city = "Iran"
#city = "Japan"
#city = "Germany"

#データを加工する
t_cases = 0

for i in range(0, len(data_r), 1):
    #if (data_r.iloc[i][1] == city): #for country/region
    if (data_r.iloc[i][0] == city):  #for province:/state  
        print(str(data_r.iloc[i][0]) + " of " + data_r.iloc[i][1])
        for day in range(4, len(data.columns), 1):            
            confirmed_r[day - 4] += data_r.iloc[i][day]
            #t_recover += data_r.iloc[i][day]
for i in range(0, len(data), 1):
    #if (data.iloc[i][1] == city): #for country/region
    if (data.iloc[i][0] == city):  #for province:/state  
        print(str(data.iloc[i][0]) + " of " + data.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed[day - 4] += data.iloc[i][day] -  confirmed_r[day - 4]            
for i in range(0, len(data_d), 1):
    #if (data_d.iloc[i][1] == city): #for country/region
    if (data_d.iloc[i][0] == city):  #for province:/state  
        print(str(data_d.iloc[i][0]) + " of " + data_d.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed_d[day - 4] += data_d.iloc[i][day] #fro drawings
            #t_deaths += data_d.iloc[i][day]
#define differencial equation of seir model
def seir_eq(v,t,beta,lp,ip,S0):
    N=N0 #int(26749*1) #58403*4
    a = -beta*v[0]*v[2]/N
    b = beta*v[0]*v[2]/N-(1/lp)*v[1] #N
    c = (1/lp)*v[1]-(1/ip)*v[2]
    d = (1/ip)*v[2]
    return [a,b,c,d]

#solve seir model
N0 = 70961 #Germany36386 #japan1482#spein46487 #korea12664 #Itary80815 #Iran35373 #58403 #korea8901 #70939. #Hubei #22116 Iran #1042 japan
N,S0,E0,I0,R0=int(N0*1.5),int(N0*1.5),int(318),int(389),0
ini_state=[S0,E0,I0,R0]
beta,lp,ip=1, 2, 7.4 #0.16049559386034234, 6.684238418998367e-32, 74.73237056193248 #Iran 0.16744196655146937, 8.355034850220927e-24, 50.50204933956246 #Hubei 0.48637152205856893, 2.4823823094760824, 27.376347305942296, 88616.56740233 #1, 2, 7.4 
t_max=len(days_from_22_Jan_20)
dt=0.01
t=np.arange(0,t_max,dt)
plt.plot(t,odeint(seir_eq,ini_state,t,args=(beta,lp,ip,N))) #0.0001,1,3
plt.legend(['Susceptible','Exposed','Infected','Recovered'])
#plt.pause(1)
#plt.close()

obs_i_2 = confirmed
obs_i_3 = confirmed_r

plt.plot(obs_i_2,"o", color="red",label = "data_2")
plt.plot(obs_i_3,"o", color="green",label = "data_3")
plt.legend()
plt.pause(1)
plt.close()

#function which estimate i from seir model func 
def estimate_i(ini_state,beta,lp,ip,N):
    v=odeint(seir_eq,ini_state,t,args=(beta,lp,ip,N))
    est=v[0:int(t_max/dt):int(1/dt)] 
    return est[:,2],est[:,3] #v0-S,v1-E,v2-I,v3-R
    
#define logscale likelihood function

def y(params):
    est_i_2,est_i_3=estimate_i(ini_state,params[0],params[1],params[2],params[3])
    return np.sum(1*(est_i_2-obs_i_2)*(est_i_2-obs_i_2)+1*(est_i_3-obs_i_3)*(est_i_3-obs_i_3))
"""
def y(params):
    est_i_2,est_i_3=estimate_i(ini_state,params[0],params[1],params[2],params[3])
    return 0*np.sum((est_i_2-obs_i_2*np.log(est_i_2)))+1*np.sum((est_i_3-obs_i_3*np.log(est_i_3))) #np.sum((est_i_2-obs_i_2*np.log(est_i_2)))
"""    
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
lns1=ax1.plot(confirmed,"o", color="red",label = "data_I")
lns2=ax1.plot(confirmed_r,"o", color="green",label = "data_R")
est_i_2,est_i_3=estimate_i(ini_state,mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3])
lns3=ax1.plot(est_i_2, label = "estimation_I")
lns0=ax1.plot(est_i_3, label = "estimation_R")

lns_ax1 = lns1+lns2 + lns3 + lns0
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)
ax1.set_ylim(0,N0)

t_max=200 #len(days_from_22_Jan_20)
t=np.arange(0,t_max,dt)

lns4=ax2.plot(confirmed,"o", color="red",label = "data")
lns5=ax2.plot(confirmed_r,"*", color="green",label = "recovered")
lns6=ax2.plot(t,odeint(seir_eq,ini_state,t,args=(mnmz.x[0],mnmz.x[1],mnmz.x[2],mnmz.x[3])))
ax2.legend(['data','data_r','Susceptible','Exposed','Infected','Recovered'], loc=1)
ax2.set_title('SEIR_b{:.2e}_ip{:.2e}_gamma{:.0f}_N{:.0f}_E0{:d}_I0{:d}_R0{:.0f}'.format(beta_const,lp,gamma_const,N,E0,I0,R0))
ax1.grid()
ax2.grid()
plt.savefig('./fig/SEIR_{}_b{:.2e}_ip{:.2e}_gamma{:.0f}_N{:.0f}_E0{:d}_I0{:d}_R0{:.0f}_.png'.format(city,beta_const,lp,gamma_const,N,E0,I0,R0)) 
plt.show()
plt.close()