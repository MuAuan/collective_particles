import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time

PARTICLE_NO = 1000-800 # 粒子数
ITERATION = 200 # 最大ループ回数　感染者が０になると止まる
MIN_X, MIN_Y = -100.0, -100.0 # 探索開始時の範囲最小値
MAX_X, MAX_Y = 100.0, 100.0 # 探索開始時の範囲最大値
recovery=30 #一定時間経過したら治癒
p=0.5*2/10 #probability of infecion
mvr = 1
rc = 400
ps0 = 1
ps1 = 1

start = time.time()

def plot_particle(sk,positions,elt,r,g,b):
    #fig, ax = plt.subplots()
    el_time = time.time()-start
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False,figsize=(16, 8))
        
    for j in range(0,PARTICLE_NO):
        x=positions[j]["x"]
        y=positions[j]["y"]
        c=positions[j]["c"]
        s = 5**2
        ax1.scatter(x, y, s, c, marker="o")
    ax1.set_xlim([MIN_X, MAX_X])
    ax1.set_ylim([MIN_Y, MAX_Y])
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title("{:.2f}:InfectionRate;{:.2f} %".format(el_time,(PARTICLE_NO-b[-1])/PARTICLE_NO*100))
    
    ind = np.arange(len(elt))  # the x locations for the groups
    width = 0.3       # the width of the bars

    ax2.set_ylim([0, PARTICLE_NO])
    ax2.set_title("{:.2f}:red_{} green_{} blue_{}".format(el_time,r[-1],g[-1],b[-1]))
    rect1 = ax2.bar(ind, b,width, color="b")
    rect2 = ax2.bar(ind+width, g, width, color="g") #, bottom=b)
    rect3 = ax2.bar(ind+2*width, r,width, color="r") #, bottom=b)
    #plt.pause(0.1)
    plt.savefig('./fig/fig{}_.png'.format(sk)) 
    plt.close()
    
# 粒子の位置更新関数
def update_position(positions,velocity):
    x0 = []
    y0 = []
    for i in range(PARTICLE_NO):
        c=positions[i]["c"]
        t_time = positions[i]["t"]  #初期値０，感染は感染時時間
        k_time = time.time()-start  #経過時間
        s = positions[i]["flag"]    #感染なし０，感染：１
        if s == 1 and c == "red":   #感染済な場合
            if k_time-t_time>recovery:  #一定時間経過したら治癒
                #print("inside",i,s,c,k_time-t_time)
                c = "blue"
                positions[i]["c"] = "green"
                positions[i]["flag"] = 1   #ただし、感染履歴ありのまま
        if c == "red":  #感染redなら位置情報取得
            x0.append(positions[i]["x"])
            y0.append(positions[i]["y"])
            #print("4",i,s,c,t_time)
    #print(x0,y0)   
    position = []
    for j in range(PARTICLE_NO):
        x=positions[j]["x"]
        y=positions[j]["y"]
        c=positions[j]["c"]
        s = positions[j]["flag"]
        t_time = positions[j]["t"]
        for k in range(len(x0)):
            if (x-x0[k])**2+(y-y0[k])**2 < rc and random.uniform(0,1)<p:
                if s ==0:
                    c = "red"
                    t_time = time.time()-start
                    s = 1
                    positions[j]["flag"]=s
                else:
                    continue
        vx = velocity[j]["x"]+mvr*random.uniform(-1, 1) #係数が粒子の運動性の大きさ
        vy = velocity[j]["y"]+mvr*random.uniform(-1, 1)
        new_x = x + vx
        new_y = y + vy
        p_color = c
        s=s

        position.append({"x": new_x, "y": new_y, "c": p_color, "t": t_time,"flag":s})
        velocity.append({"x": vx, "y": vy})

    return position, velocity, x0

def count_brg(position):
    r=0
    g=0
    b=0
    for j in range(len(position)):
        if position[j]["c"] == "red":
            r += 1
        elif position[j]["c"] == "green":
            g += 1
        else:
            b += 1
    return r,g,b        

def main():
    # 時間計測開始
    #start = time.time()
    xy_min, xy_max = -32, 32
    # 各粒子の初期位置, 速度, personal best, global best 及びsearch space設定
    position = []
    velocity = []  #速度は使えるように拡張

    # 初期位置, 初期速度
    for s in range(0,ps0):
        #position.append({"x": random.uniform(MIN_X, MAX_X), "y": random.uniform(MIN_Y, MAX_Y), "c": "red", "t":0, "flag":1})
        position.append({"x": 0*random.uniform(MIN_X, MAX_X), "y": 0*random.uniform(MIN_Y, MAX_Y), "c": "red", "t":0, "flag":1}) #真ん中（0，0）に初期感染者を1人置く
        velocity.append({"x": 0, "y": 0}) #感染者の初速度0としている
    for s in range(ps0,ps1):
        #position.append({"x": random.uniform(MIN_X, MAX_X), "y": random.uniform(MIN_Y, MAX_Y), "c": "red", "t":0, "flag":1})
        position.append({"x": random.uniform(MIN_X, MAX_X), "y": random.uniform(MIN_Y, MAX_Y), "c": "green", "t":0, "flag":1}) #真ん中（0，0）に初期感染者を1人置く
        velocity.append({"x": 0, "y": 0}) #感染者の初速度0としている        
    for s in range(ps1,PARTICLE_NO):
        position.append({"x": random.uniform(MIN_X, MAX_X), "y": random.uniform(MIN_Y, MAX_Y), "c": "blue", "t": 0, "flag":0})
        velocity.append({"x": 0, "y": 0})
    """    
    for i in range(0,10): #x方向に集落を10個並べる
        for j in range(0,10): #xy方向に集落をメッシュ10で並べる
            for k in range(0,10): #１集落辺り10個の感受性保持者を分布
                s=k+j*10+i*100;
                position.append({"x": 10+(-100+i*20)+random.uniform(MIN_X/100, MAX_X/100), "y":10+(-100+j*20)+ random.uniform(MIN_Y/100, MAX_Y/100), "c": "blue", "t": 0, "flag":0})
                velocity.append({"x": 0, "y": 0})
    """            
    print(len(position))
    sk = 0
    red=[]
    green=[]
    blue=[]
    elapsed_time = []
    while sk < ITERATION:
        position,velocity, x0 = update_position(position,velocity) ######
        r,g,b = count_brg(position)
        red.append(r)
        green.append(g)
        blue.append(b)
        el_time=time.time()-start
        elapsed_time.append(el_time)
        #print("{:.2f}:red_{} green_{} blue_{}".format(el_time,r,g,b))
        plot_particle(sk,position,elapsed_time,red,green,blue)
        if x0==[]:
            break
        sk += 1

    # 時間計測終了
    process_time = time.time() - start
    print("time:", process_time)

    
if __name__ == '__main__':
    main()    
    