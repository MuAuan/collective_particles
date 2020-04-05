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
day_confirmed = [0] * (len(data.columns) - 4)
confirmed_r = [0] * (len(data_r.columns) - 4)
day_confirmed_r = [0] * (len(data.columns) - 4)
confirmed_d = [0] * (len(data_d.columns) - 4)
diff_confirmed = [0] * (len(data.columns) - 4)
days_from_22_Jan_20 = np.arange(0, len(data.columns) - 4, 1)
days_from_22_Jan_20_ = np.arange(0, len(data.columns) - 4, 1)
beta_ = [0] * (len(data_r.columns) - 4)
gamma_ = [0] * (len(data_d.columns) - 4)


city = "Japan"

skd=5 #2 #1 #4 #3 #2 #slopes average factor
#データを加工する
t_cases = 0
t_recover = 0
t_deaths = 0
for i in range(0, len(data_r), 1):
    if (data_r.iloc[i][1] == city): #for country/region
    #if (data_r.iloc[i][0] == city):  #for province:/state  
        print(str(data_r.iloc[i][0]) + " of " + data_r.iloc[i][1])
        for day in range(4, len(data.columns), 1):            
            confirmed_r[day - 4] += data_r.iloc[i][day]
            if day < 4+skd:
                day_confirmed_r[day-4] += data_r.iloc[i][day]
            else:
                day_confirmed_r[day-4] += (data_r.iloc[i][day] - data_r.iloc[i][day-skd])/(skd)
        t_recover += data_r.iloc[i][day]        
for i in range(0, len(data_d), 1):
    if (data_d.iloc[i][1] == city): #for country/region
    #if (data_d.iloc[i][0] == city):  #for province:/state  
        print(str(data_d.iloc[i][0]) + " of " + data_d.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed_d[day - 4] += data_d.iloc[i][day] #fro drawings
        t_deaths += data_d.iloc[i][day]        
for i in range(0, len(data), 1):
    if (data.iloc[i][1] == city): #for country/region
    #if (data.iloc[i][0] == city):  #for province:/state  
        print(str(data.iloc[i][0]) + " of " + data.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed[day - 4] += data.iloc[i][day] -  confirmed_r[day - 4] -confirmed_d[day-4]
            diff_confirmed[day - 4] += confirmed[day-4] /  (confirmed_r[day - 4]+confirmed_d[day-4])
            if day == 4:
                day_confirmed[day-4] += data.iloc[i][day]
            else:
                day_confirmed[day-4] += data.iloc[i][day] - data.iloc[i][day-1]

tl_confirmed = 0
dlog_confirmed = [0] * (len(data.columns) - 4)
dlog_confirmed[0]=np.log(confirmed[0])
dlog_confirmed[1]=np.log(confirmed[1])-np.log(confirmed[0])
ratio_confirmed = [0] * (len(data.columns) - 4)
ratio_confirmed[0]=np.log(confirmed[0])
ratio_confirmed[1]=(confirmed[1]-confirmed[0])/(confirmed[0])
ratio_confirmed[2]=(confirmed[2]-confirmed[0])/(confirmed[0])/2

for i in range(skd, len(confirmed), 1):        
    if confirmed[i] > 0:    
        gamma_[i]=day_confirmed_r[i]/confirmed[i]
    else:
        continue
tl_confirmed = confirmed[len(confirmed)-1] + confirmed_r[len(confirmed)-1] + confirmed_d[len(confirmed)-1]
t_cases = tl_confirmed

t_max=len(confirmed)
dt=1
t=np.arange(0,t_max,dt)
t1=t

obs_i = confirmed_r[:]
#function which estimate i from seir model func 
def estimate_i(ini_state,r0,alpha):
    est = r0*np.exp(alpha*(t))
    return est

def y(params):
    est_i=estimate_i(ini_state,params[0],params[1])
    return np.sum((est_i-obs_i)*(est_i-obs_i))
r0=1
alpha = 1
ini_state=[5.70579672, 0.00755685]
#optimize logscale likelihood function
mnmz=minimize(y,ini_state,method="nelder-mead")
print(mnmz)
r0,alpha = mnmz.x[0],mnmz.x[1]
est=estimate_i(ini_state,r0,alpha)

t=np.arange(63,t_max,dt)
t2=t
obs_i = confirmed[63:]
r0_=1
alpha_ = 1
ini_state=[5.70579672, 0.00755685]
#optimize logscale likelihood function
mnmz=minimize(y,ini_state,method="nelder-mead")
print(mnmz)
r0_,alpha_ = mnmz.x[0],mnmz.x[1]
#est_confirmed=estimate_i(ini_state,r0_,alpha_)
#t=np.arange(63,100,dt)
t3=t
est_confirmed=estimate_i(ini_state,r0_,alpha_)

diff_est=[0] * (len(data.columns) - 4)
gamma_est=[0] * (len(data.columns) - 4)
R_est = [0] * (len(data_d.columns) - 4)
R_0 = [0] * (len(data_d.columns) - 4)
for i in range(1,t_max):
    diff_est[i]=est[i]-est[i-1]
for i in range(0, len(confirmed), 1):        
    if confirmed[i] > 0:    
        gamma_est[i]=diff_est[i]/confirmed[i]
        R_est[i]= 1+day_confirmed[i]/diff_est[i] # diff_est=gamma*confirmed
        R_0[i]= R_est[i]/(1-gamma_est[i]*R_est[i]*confirmed[i]*i/t_cases)
    else:
        continue
    
#matplotlib描画
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
#ax3 = ax1.twinx()
ax4 = ax2.twinx()

lns1=ax1.semilogy(days_from_22_Jan_20, confirmed, "o-", color="red",label = "cases")
lns8=ax1.semilogy(t3, est_confirmed, "-", color="black",label = "cases_r0_={:.2f}alpha_={:.2e}".format(r0_,alpha_))
lns2=ax1.semilogy(days_from_22_Jan_20, confirmed_r, "*-", color="green",label = "recovered+deaths")
#lns4=ax2.plot(days_from_22_Jan_20_, dlog_confirmed, "o-", color="blue",label = "dlog_confirmed")
#lns3=ax4.plot(days_from_22_Jan_20_, gamma_, "o-", color="black", zorder=1,label = "gamma")
lns3=ax4.plot(days_from_22_Jan_20_, gamma_est, "o-", color="black", zorder=1,label = "gamma_est")
#lns4=ax2.bar(days_from_22_Jan_20_, day_confirmed, zorder=2, label = "day_confirmed")
lns4=ax2.plot(days_from_22_Jan_20_, R_est, "o-", color="blue",label = "R_est")
lns5=ax1.semilogy(days_from_22_Jan_20_, diff_confirmed, ".-", color="black",label = "I/(R+D)")
lns6=ax1.semilogy(t1, est,"-", color="black", zorder=1, label = "est_r0={:.2f}alpha={:.2e}".format(r0,alpha))
lns7=ax2.plot(t1, diff_est,"-", color="black", zorder=1, label = "diff_est_r0={:.2f}alpha={:.2e}".format(r0,alpha))
lns9=ax2.bar(days_from_22_Jan_20_, day_confirmed_r, label = "day_confirmed_r")
#lns10=ax2.plot(days_from_22_Jan_20_, R_0, "o-", color="red",label = "R_0")

lns_ax1 = lns1 +lns2 +lns5 + lns6 +lns8
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)

lns_ax2 = lns3 #+lns9
labs_ax2 = [l.get_label() for l in lns_ax2]
ax4.legend(lns_ax2, labs_ax2, loc=0)
ax2.legend(loc=2)

ax1.set_title(city +" ; {} cases, {} recovered, {} deaths".format(t_cases,t_recover,t_deaths))
ax1.set_xlabel("days from 22, Jan, 2020")
ax1.set_ylabel("casas, recovered ")
#ax2.set_ylabel("dlog_confirmed")
ax4.set_ylabel("gamma")
ax2.set_ylabel("day_confirmed_r, R")
ax4.set_ylim(0,0.04)
ax2.set_ylim(0,40)
#ax2.set_yscale('log')
#ax4.set_yscale('log')

#ax3.set_ylabel("deaths ")
#ax4.set_ylabel("deaths_rate %")
#ax4.set_ylim(-0.5,0.5)
ax1.grid()
ax2.grid()

plt.pause(1)
plt.savefig('./fig/removed_{}_gamma_R_{}.png'.format(city,skd)) 
plt.close() 

t=np.arange(63,100,dt)
t3=t
est_confirmed=estimate_i(ini_state,r0_,alpha_)

#matplotlib描画
fig, ax3 = plt.subplots(1,1,figsize=(1.6180 * 4, 4*1))

lns1=ax3.semilogy(days_from_22_Jan_20, confirmed, "o-", color="red",label = "cases")
lns8=ax3.semilogy(t3, est_confirmed, "-", color="black",label = "cases_r0_={:.2f}alpha_={:.2e}".format(r0_,alpha_))
lns2=ax3.semilogy(days_from_22_Jan_20, confirmed_r, "*-", color="green",label = "recovered+deaths")
lns5=ax3.semilogy(days_from_22_Jan_20_, diff_confirmed, ".-", color="black",label = "I/(R+D)")
lns6=ax3.semilogy(t1, est,"-", color="black", zorder=1, label = "est_r0={:.2f}alpha={:.2e}".format(r0,alpha))

lns_ax1 = lns1 +lns2 +lns5 + lns6 +lns8
labs_ax1 = [l.get_label() for l in lns_ax1]
ax3.legend(lns_ax1, labs_ax1, loc=0)

ax3.set_title(city +" ; {} cases, {} recovered, {} deaths".format(t_cases,t_recover,t_deaths))
ax3.set_xlabel("days from 22, Jan, 2020")
ax3.set_ylabel("casas, recovered ")
ax3.grid()

plt.pause(1)
plt.savefig('./fig/exterpolate_{}_gamma_R_{}.png'.format(city,skd)) 
plt.close() 