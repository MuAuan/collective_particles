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
#city = "Hubei"
#city = "Korea, South"
#city = "Italy"
#city = "Spain"
#city = "Iran"
#city = "Japan"
#city = "Germany"
#city = "Bahrain"
#city = "Switzerland"
city = "United Kingdom"

#データを加工する
t_cases = 0

for i in range(0, len(data_r), 1):
    if (data_r.iloc[i][1] == city): #for country/region
    #if (data_r.iloc[i][0] == city):  #for province:/state  
        print(str(data_r.iloc[i][0]) + " of " + data_r.iloc[i][1])
        for day in range(4, len(data.columns), 1):            
            confirmed_r[day - 4] += data_r.iloc[i][day]
            #t_recover += data_r.iloc[i][day]
for i in range(0, len(data), 1):
    if (data.iloc[i][1] == city): #for country/region
    #if (data.iloc[i][0] == city):  #for province:/state  
        print(str(data.iloc[i][0]) + " of " + data.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed[day - 4] += data.iloc[i][day] -  confirmed_r[day - 4]            
for i in range(0, len(data_d), 1):
    if (data_d.iloc[i][1] == city): #for country/region
    #if (data_d.iloc[i][0] == city):  #for province:/state  
        print(str(data_d.iloc[i][0]) + " of " + data_d.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed_d[day - 4] += data_d.iloc[i][day] #fro drawings
            #t_deaths += data_d.iloc[i][day]
            
#define differencial equation of sir model
def sir_eq(v,t,beta,gamma):
    a = -beta*v[0]*v[1]
    b = beta*v[0]*v[1] - gamma*v[1]
    c = gamma*v[1]
    return [a,b,c]

def estimate_i(ini_state,beta,gamma):
    v=odeint(sir_eq,ini_state,t,args=(beta,gamma))
    est=v[0:int(t_max/dt):int(1/dt)] 
    return est[:,0],est[:,1],est[:,2] #v0-S,v1-E,v2-I,v3-R

#define logscale likelihood function
def y(params):
    est_i_0,est_i_1,est_i_2=estimate_i(ini_state,params[0],params[1])
    return np.sum(f1*(est_i_1-obs_i_2)*(est_i_1-obs_i_2)+f2*(est_i_2-obs_i_3)*(est_i_2-obs_i_3))

#solve seir model
N0 = 14600 #13159 #51213 #70857 #70975 #470 #34710 #9471 #1517
f1,f2 = 1,0 #data & data_r fitting factor
N,S0,I0,R0=int(N0*50),int(N0*50),int(1000),int(0)
ini_state=[S0,I0,R0]
beta,gamma = 1.6e-6, 1e-3 #2.3400696465302977e-06, 0.053333851016341804 #4.5292086907047916e-05, 0.015805239256950265

t_max=len(days_from_22_Jan_20)
dt=0.01
t=np.arange(0,t_max,dt)
plt.plot(t,odeint(sir_eq,ini_state,t,args=(beta,gamma))) #0.0001,1,3
plt.legend(['Susceptible','Infected','Recovered'])

obs_i_2 = confirmed
obs_i_3 = confirmed_r

plt.plot(obs_i_2,"o", color="red",label = "data_2")
plt.plot(obs_i_3,"o", color="green",label = "data_3")
plt.legend()
plt.pause(1)
plt.close()

#optimize logscale likelihood function
mnmz=minimize(y,[beta,gamma],method="nelder-mead")
print(mnmz)
#R0
#N_total = S_0+I_0+R_0
#R0 = N_total*beta_const *(1/gamma_const)
beta_const,gamma = mnmz.x[0],mnmz.x[1] #感染率、除去率（回復率）
print(beta_const,gamma)
r0 = N*beta_const*(1/gamma)
print(r0)

#plot reult with observed data
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
ax3 = ax1.twinx()
ax4 = ax2.twinx()

lns1=ax1.plot(confirmed,"o", color="red",label = "data_I")
lns2=ax1.plot(confirmed_r,"o", color="green",label = "data_R")
est_i_0,est_i_1,est_i_2=estimate_i(ini_state,mnmz.x[0],mnmz.x[1])
lns3=ax1.plot(est_i_1, label = "estimation_I")
lns0=ax1.plot(est_i_2, label = "estimation_R")
lns_1=ax3.plot((S0-est_i_1-est_i_2)*r0/N,".", color="black", label = "effective_R0")
ax3.set_ylim(0,r0+1)

lns_ax1 = lns1+lns2 + lns3 + lns0 +lns_1
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)
ax1.set_ylim(0,N0)
ax1.set_title('SIR_{} f1_{:,d} f2_{:,d};N={:.0f} S0={:.0f} I0={:.0f} R0={:.0f} R={:.2f}'.format(city,f1,f2,N,S0,I0,R0,(S0-est_i_1[-1]-est_i_2[-1])*r0/N))
ax1.set_ylabel("Susceptible, Infected, recovered ")
ax3.set_ylabel("effective_R0")

t_max=200 #len(days_from_22_Jan_20)
t=np.arange(0,t_max,dt)

lns4=ax2.plot(confirmed,"o", color="red",label = "data")
lns5=ax2.plot(confirmed_r,"o", color="green",label = "data_r")
est_i_0,est_i_1,est_i_2=estimate_i(ini_state,mnmz.x[0],mnmz.x[1])
lns6=ax2.plot(est_i_0, label = "estimation_S")
lns7=ax2.plot(est_i_1, label = "estimation_I")
lns8=ax2.plot(est_i_2, label = "estimation_R")
lns_6=ax4.plot((S0-est_i_1-est_i_2)*r0/N,".", color="black", label = "effective_R0")

lns_ax2 = lns4+lns5 + lns6 + lns7+ lns8 +lns_6
labs_ax2 = [l.get_label() for l in lns_ax2]
ax2.legend(lns_ax2, labs_ax2, loc=0)
ax4.set_ylim(0,r0+1)
ax2.set_ylim(0,S0)
ax2.set_title('SIR_{};b={:.2e} g={:.2e} r0={:.2f}'.format(city,beta_const,gamma,r0))
ax2.set_ylabel("Susceptible, Infected, recovered ")
ax4.set_ylabel("effective_R0")

ax1.grid()
ax2.grid()
plt.savefig('./fig/SIR_{}f1_{:,d}f2_{:,d};b_{:.2e}g_{:.2e}r0_{:.2f}N_{:.0f}S0_{:.0f}I0_{:.0f}R0_{:.0f}.png'.format(city,f1,f2,beta_const,gamma,r0,N,S0,I0,R0)) 
plt.show()
plt.close()
