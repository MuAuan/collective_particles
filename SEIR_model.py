import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define differential equation of SEIR model

'''
dS/dt = -beta * S * I / N
dE/dt = beta* S * I / N - epsilon * E
dI/dt = epsilon * E - gamma * I
dR/dt = gamma * I

[v[0], v[1], v[2], v[3]]=[S, E, I, R]

dv[0]/dt = -beta * v[0] * v[2] / N
dv[1]/dt = beta * v[0] * v[2] / N - epsilon * v[1]
dv[2]/dt = epsilon * v[1] - gamma * v[2]
dv[3]/dt = gamma * v[2]

'''


def SEIR_EQ(v, t, beta, epsilon, gamma, N ):
    a=-beta * v[0] * v[2] / N
    b=beta * v[0] * v[2] / N - epsilon * v[1]
    c=epsilon * v[1] - gamma * v[2]
    d=gamma * v[2]
    return [a,b,c,d]
    #return [-beta * v[0] * v[2] / N ,beta * v[0] * v[2] / N - epsilon * v[1],
            #epsilon * v[1] - gamma * v[2],gamma * v[2]]
            
  # parameters
t_max = 60 #100 #days
dt = 0.01

# initial_state
S_0 = 80000 #14000000 #99
E_0 = 618 #318 #1
I_0 = 180 #389 #389 #0
R_0 = 0
N_pop = S_0 + E_0 + I_0 + R_0
ini_state = [S_0, E_0, I_0, R_0]  # [S[0],E,[0], I[0], R[0]]


#感染率
beta_const = 0.75 #1 # 0.18  #1 #感染率

#暴露後に感染症を得る率
latency_period = 5 #2 #2 #days
epsilon_const = 1/latency_period

#回復率や隔離率
infectious_period = 30 #14 #7.4 #7.4 #days
gamma_const = 1/infectious_period

# numerical integration
times = np.arange(0, t_max, dt)
args = (beta_const, epsilon_const, gamma_const, N_pop)

# Numerical Solution using scipy.integrate
# Solver SEIR model
result = odeint(SEIR_EQ, ini_state, times, args)
S = []
E = []
I = []
R = []
R_rate = []
il_p=0
rl_p=0
r_rate=0
for i in range(len(result)):
    sl = result[i][0]
    el = result[i][1]
    il = result[i][2]
    rl = result[i][3]
    il_p = il
    rl_p = rl
    r_rate = rl_p*100/(il_p+rl_p)
    S.append(sl)
    E.append(el)
    I.append(il)
    R.append(rl)
    R_rate.append(r_rate)
    
print(R_rate[1000:1010])
print( R[1000:1010])
print(I[1000:1010])
# plot
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
#ax1.plot(times, S, label='Susceptible')
ax1.plot(times, R, label ='Recovered')
#ax2.plot(times, E, label='Exposed')
ax1.plot(times, I, label ='Infectious')
ax2.plot(times, R_rate, label ='recovered')
ax1.legend() #['Susceptible', 'Removed']) #, 'Exposed']) #, 'Infectious', 'Removed'])
ax2.legend()
ax1.set_title("SEIR model  COVID-19")
ax1.set_xlabel('time(days)')
ax1.set_ylabel('population')
ax2.set_ylabel('recovered_rate %')
ax1.grid()
ax2.grid()

plt.show()