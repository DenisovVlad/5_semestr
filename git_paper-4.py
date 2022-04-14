
import numpy as np
from scipy.integrate import RK45
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random

##################
##### RUNGE-KUTTA
##################
# The following utility returns the function dP_dt given two vectors of coefficients $a$ and $b$ of equal size
def get_dP_dt(a, b):
    assert isinstance(a, np.ndarray)
    assert isinstance(b, np.ndarray)
    assert len(a.shape) == len(b.shape) == 1
    assert a.shape[0] > 0
    assert a.shape[0] == b.shape[0]

    dim = a.shape[0]
    def dP_dt(t,P):
        #print(P)
        J = a * P[0] * P
        J[:-1] -= b[1:] * P[1:]
        #J = np.zeros(n)
        J[-1] = 0

        c = np.empty(dim)
        c[0] = -J[0] - J.sum()
        c[1:] = J[:-1] - J[1:]
        return c
    
    return dP_dt

n = 10
# The parameters $a$ and $b$ are the vectors containing all coefficients $a_i$ and $b_i$. They must contain $n$ elements.
# In this moment all coefficients are $1$.

#k_1 = random.randint(1,1001)
#k_2 = random.randint(1,1001)
k_1=1
k_2=1
dP_dt = get_dP_dt(
    #a=np.ones(n),
    #b=np.ones(n)
    a= np.full(n, k_1),
    b= np.full(n, k_2),
)

# Time and Starting point coordinates


time= 1000
n_i=time-1
#ts = np.linspace(0, 1500, time)
ts = 50
#P0 = np.array([2.0, 0.0])

P0 = np.zeros(n)   # inizialize an array with n elements
P0[0] = n
#print (P0)
# Resolution
solution = RK45(dP_dt, t0 = 0, y0 = P0, t_bound = ts)

# collect data
t_values = []
y_values = []

for i in range(n):
    y_values.append([])

for i in range(ts):
    # get solution step state
    solution.step()
    t_values.append(solution.t)
    
    y_values[0].append(solution.y[0])
    y_values[1].append(solution.y[1])
    y_values[2].append(solution.y[2])
    y_values[3].append(solution.y[3])
    y_values[4].append(solution.y[4])
    y_values[5].append(solution.y[5])
    y_values[6].append(solution.y[6])
    y_values[7].append(solution.y[7])
    y_values[8].append(solution.y[8])
    y_values[9].append(solution.y[9])
    
    
    # break loop after modeling is finished
    if solution.status == 'finished':
        break
#print(t_values)
#print(y_values)

mean_size = []
for j in range(ts):
    mean_size.append([])
    for i in range(10):
        mean_size[j].append(0)
res = []
for j in range(ts):
    res.append(0)
for i in range(ts):
    for j in range(1,5):
        res[i] += y_values[j][i]*(j)
    res[i] = res[i]/4



# # PLOT
plt.figure(figsize=(10,5))
for col in range(1):
    plt.plot(t_values, y_values[col], label="c{}".format(col+1))
plt.xlabel("Time")
plt.ylabel("Concentration of single paticles")
plt.legend()
plt.xlim(-0.2, ts/2)
plt.ylim(0., 5)
#plt.xlim(0.3, ts/2)
#plt.ylim(0.02, 1.4)
plt.show()

# # PLOT
plt.figure(figsize=(10,5))
for col in range(1,5):
    plt.plot(t_values, y_values[col], label="c{}".format(col+1))
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.xlim(-0.2, ts/2)
plt.ylim(0., 5)
#plt.xlim(0.3, ts/2)
#plt.ylim(0.02, 1.4)
plt.show()

# # PLOT
plt.figure(figsize=(10,5))
plt.plot(t_values, res)
plt.xlabel("Time")
plt.ylabel("Mean size")
plt.legend()
plt.xlim(0.07, 5)
#plt.ylim(1, 1.5)
plt.show()




'''
#######################
## ODEINT
def get_dP_dt1(a, b):
    assert isinstance(a, np.ndarray)
    assert isinstance(b, np.ndarray)
    assert len(a.shape) == len(b.shape) == 1
    assert a.shape[0] > 0
    assert a.shape[0] == b.shape[0]
    
    dim = a.shape[0]
    def dP_dt1(P, t):
        J = a * P[0] * P
        J[:-1] -= b[1:] * P[1:]
        J[-1] = 0
        
        c = np.empty(dim)
        c[0] = -J[0] - J.sum()
        c[1:] = J[:-1] - J[1:]
        return c
    
    return dP_dt1

n = 5


# The parameters $a$ and $b$ are the vectors containing all coefficients $a_i$ and $b_i$. They must contain $n$ elements.
# In this moment all coefficients are $1$.

#k_1 = random.randint(1,1001)
#k_2 = random.randint(1,1001)
k_1=1
k_2=1
dP_dt1 = get_dP_dt1(
                  #a=np.ones(n),
                  #b=np.ones(n)
                  a= np.full(n, k_1),
                  b= np.full(n, k_2),
                  )


# Time and Starting point coordinates

time= 1000
n_i=time-1
ts = np.linspace(0, 500, time)
#P0 = np.array([2.0, 0.0])

P0 = np.zeros(n)   # inizialize an array with n elements
P0[0] = n


Ps = odeint(dP_dt1, P0, ts)

c1 = Ps[:,0]
#print c1[n_i]


# # PLOT
plt.figure(figsize=(15,5))
for col in range(n):
    plt.plot(ts, Ps[:,col], label="c{}".format(col+1))
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
'''
