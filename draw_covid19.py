# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#pandasでCSVデータ読む。
data = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv')
data_r = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv')
data_d = pd.read_csv('COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv')

confirmed = [0] * (len(data.columns) - 4)
confirmed_r = [0] * (len(data_r.columns) - 4)
confirmed_d = [0] * (len(data_d.columns) - 4)
recovered_rate = [0] * (len(data_r.columns) - 4)
deaths_rate = [0] * (len(data_d.columns) - 4)
days_from_22_Jan_20 = np.arange(0, len(data.columns) - 4, 1)

city = "Hubei"
city = "Korea, South"
city = "Iran"
city = "Italy"
city = "Spain"
city = "Iraq"
city = "Japan"
city = "Singapore"
city = "Germany"
city = "China"
city = "US"
city = "France"
city = "United Kingdom"
city = "Switzerland"
city = "Indonesia"
#city = "Guangdong"
#city = "Zhejiang"
#city = "New York"

#データを加工する
t_cases = 0
t_recover = 0
t_deaths = 0
for i in range(0, len(data), 1):
    if (data.iloc[i][1] == city): #for country/region
    #if (data.iloc[i][0] == city):  #for province:/state  
        print(str(data.iloc[i][0]) + " of " + data.iloc[i][1])
        for day in range(4, len(data.columns), 1):
            confirmed[day - 4] += data.iloc[i][day] -  data_r.iloc[i][day] #+=
            confirmed_r[day - 4] += data_r.iloc[i][day]
            confirmed_d[day - 4] += data_d.iloc[i][day]*1 #fro drawings
 
        #t_cases += data.iloc[i][day] 
        t_recover += data_r.iloc[i][day]        
        t_deaths += data_d.iloc[i][day]
tl_confirmed = 0        
for i in range(0, len(confirmed), 1):
    tl_confirmed = confirmed[i] + confirmed_r[i] + confirmed_d[i]
    #print(tl_confirmed)
    if tl_confirmed > 0:
        recovered_rate[i]=float(confirmed_r[i]*100)/float(tl_confirmed)
        deaths_rate[i]=float(confirmed_d[i]*100)/float(tl_confirmed)
    else:
        continue
t_cases = tl_confirmed       
        
#print(days_from_22_Jan_20)
#print(confirmed)
#print(confirmed_r)
#print(recovered_rate)
#print(deaths_rate)

#matplotlib描画
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
ax3 = ax1.twinx()
ax4 = ax2.twinx()

lns1=ax1.plot(days_from_22_Jan_20, confirmed, "o-", color="red",label = "cases")
lns2=ax1.plot(days_from_22_Jan_20, confirmed_r, "*-", color="green",label = "recovered")
lns3=ax3.plot(days_from_22_Jan_20, confirmed_d, "D-", color="black", label = "deaths")
lns4=ax2.plot(days_from_22_Jan_20, recovered_rate, "*-", color="green",label = "recovered")
lns5=ax4.plot(days_from_22_Jan_20, deaths_rate, "D-", color="black", label = "deaths")

lns_ax1 = lns1+lns2+lns3
labs_ax1 = [l.get_label() for l in lns_ax1]
ax1.legend(lns_ax1, labs_ax1, loc=0)

lns_ax2 = lns4+lns5
labs_ax2 = [l.get_label() for l in lns_ax2]
ax2.legend(lns_ax2, labs_ax2, loc=0)

ax1.set_title(city +" ; {} cases, {} recovered, {} deaths".format(t_cases,t_recover,t_deaths))
ax1.set_xlabel("days from 22, Jan, 2020")
ax1.set_ylabel("casas, recovered ")
ax2.set_ylabel("recovered_rate %")
ax2.set_ylim(0,100)

ax3.set_ylabel("deaths ")
ax4.set_ylabel("deaths_rate %")
ax4.set_ylim(0,10)
ax1.grid()
ax2.grid()

plt.pause(1)
plt.savefig('./fig/fig_{}_.png'.format(city)) 
plt.close()