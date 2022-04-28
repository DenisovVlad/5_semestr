
import numpy as np
from scipy.integrate import RK45
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random


'''
def get_dP_dt(a, b):
    dim = a.shape[0]
    def dP_dt(P,t):
        #print(P)
        J = a * P[0] * P
        J[:-1] -= b[1:] * P[1:]
        J[0] = P[0] * P[0]
        J[0] -= P[1]
        #J = np.zeros(n)
        J[-1] = 0
        c = np.empty(dim)
        c[0] = -J[0] - J.sum()
        c[1:] = J[:-1] - J[1:]
        return c
    return dP_dt
'''
def get_dP_dt(a, b):

    dim = a.shape[0]
    def dP_dt(P,t):
        c = np.zeros(dim)
        a[0] = 2
        b[0] = 2
        #a[1] = 1
        #b[1] = 1

        c[0] = -a[0]*P[0]*(P[0]+P[1]+P[2])+a[1]*(P[1]+P[2]+P[3]+P[1])-b[0]*P[0]*P[1]+b[1]*P[2]
        c[1] = -b[0]*P[1]*(P[0]+P[1])+b[1]*(P[1]+P[2]+P[3]+P[3])+a[0]*P[0]*(P[0]-P[1])-a[1]*(P[1]-P[2])
        
        for j in range(2,n-2):
            c[j] = a[0]*(P[j-1]*P[0]-P[j]*P[0])+a[1]*(P[j+1]-P[j])+b[0]*(P[j-2]*P[1]-P[j]*P[1])+b[1]*(P[j+2]-P[j])
        #c[2] = a[0]*(P[1]*P[0]-P[2]*P[0])+a[1]*(P[3]-P[2])+b[0]*(P[0]*P[1]-P[2]*P[1])+b[1]*(P[4]-P[2])
        return c
    
    return dP_dt

n = 10



#k_1 = random.randint(1,1001)
#k_2 = random.randint(1,1001)
#k_1=1
#k_2=1
a= np.full(n, 0.0)
b= np.full(n, 0.0)
a[0] = 0.053
a[1] = 0.0000061
b[0] = 10.2
b[1] = 0.56

dP_dt = get_dP_dt(
    #a=np.ones(n),
    #b=np.ones(n)
    a,
    b,
)




time= 1000
n_i=time-1
#ts = np.linspace(0, 1500, time)
ts = 25
#P0 = np.array([2.0, 0.0])

P0 = np.zeros(n)
P0[0] = n
#print (P0)
# Resolution
#solution = RK45(dP_dt, t0 = 0, y0 = P0, t_bound = ts)


def RK4(f, u0, t0, tf , nn):
    t = np.linspace(t0, tf, nn+1)
    u = np.array((nn+1)*[u0])
    h = t[1]-t[0]
    for i in range(nn):
        k1 = h * f(u[i], t[i])
        k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
        k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
        k4 = h * f(u[i] + k3, t[i] + h)
        u[i+1] = u[i] + (k1 + 2*(k2 + k3 ) + k4) / 6
    return t,u
ts,Ps = RK4(dP_dt, P0, 0,  ts, time)
print(Ps)


# # PLOT
#ts1 = np.linspace(0, 500, time)
plt.figure(figsize=(15,5))
for col in range(n-2):
    plt.plot(ts, Ps[:,col], label="c{}".format(col+1))
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.ylim(0,n/2)
plt.legend()
plt.show()


