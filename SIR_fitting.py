import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define differential equation of SEIR model

'''
dS/dt = -beta * S * I / N
dI/dt = -beta * S * I / N - gamma * I
dR/dt = gamma * I

[v[0], v[1], v[2]]=[S, I, R]

dv[0]/dt = -beta * v[0] * v[1] / N
dv[1]/dt = beta * v[0] * v[1] / N - gamma * v[1]
dv[2]/dt = gamma * v[1]
'''

def SIR_EQ(v, t, beta, gamma, N ):
    a=-beta * v[0] * v[1] / 1   #N
    b=beta * v[0] * v[1] / 1 - gamma * v[1] #N
    c=gamma * v[1]
    return [a,b,c]
            
  # parameters
t_max = 20 #60 #100 #days
dt = 0.01 #0.01

def SIR_calc(b,ip,s0,I0,R0):
    # initial_state
    S_0 = s0 #14000000 #99
    I_0 = I0 #389 #389 #0
    R_0 = R0
    N_pop = S_0 + I_0 + R_0
    ini_state = [S_0, I_0, R_0] 
    R0 = N_pop*b *(ip)
    print(R0)

    #感染率
    beta_const = b #0.75 #1 # 0.18  #1 #感染率

    #回復率や隔離率
    infectious_period = ip #30 #14 #7.4 #7.4 #days
    gamma_const = 1/infectious_period

    # numerical integration
    times = np.arange(0, t_max, dt)
    args = (beta_const, gamma_const, N_pop)

    # Numerical Solution using scipy.integrate
    # Solver SEIR model
    result = odeint(SIR_EQ, ini_state, times, args)
    S = []
    I = []
    R = []
    R_rate = []
    rl_p=0
    r_rate=0
    for i in range(len(result)):
        sl = result[i][0]
        il = result[i][1]
        rl = result[i][2]
        il_p = il
        rl_p = rl
        r_rate = rl_p*100/(il_p+rl_p)
        S.append(sl)
        I.append(il)
        R.append(rl)
        R_rate.append(r_rate)
    data_influ=[3,8,28,75,221,291,255,235,190,125,70,28,12,5]
    data_day = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    
    # plot
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(1.6180 * 4, 4*2))
    #ax1.plot(times, S, label='Susceptible')
    ax1.plot(times, R, label ='Recovered')
    #ax2.plot(times, E, label='Exposed')
    ax1.plot(times, I, label ='Infectious')
    ax1.plot(data_day,data_influ,"o", color="red",label = 'data')
    ax2.plot(times, R_rate, label ='recovered')
    ax1.legend() #['Susceptible', 'Removed']) #, 'Exposed']) #, 'Infectious', 'Removed'])
    ax2.legend()
    ax1.set_title("SIR model  COVID-19")
    ax1.set_xlabel('time(days)')
    ax1.set_ylabel('population')
    ax2.set_ylabel('recovered_rate %')
    ax1.grid()
    ax2.grid()

    plt.pause(1)
    plt.savefig('./fig/SIR_b{}_ip{}_s0{}_I0{}_R0{}_.png'.format(b,ip,s0,I0,R0)) 
    plt.close()

def main():
    ip=2 #gamma = 0.5
    s0=762 #762
    I0=1
    R0=0
    b=0.0021 #0.0022 #0.0026*(s0+I0+R0)
    SIR_calc(b,ip,s0,I0,R0)
    for j in range(18,28,1):
        ip=float(j/10)      
        SIR_calc(b,ip,s0,I0,R0)
    
if __name__ == '__main__':
    main()        
    