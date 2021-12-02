
import matplotlib.pyplot as plt
from matplotlib import rc

import math as m
import numpy as np
import random

# Rate function of coagulation
# Скоростная функция коагуляции
def a(i):

	return i
	
# Rate function of fragmentation
# Скоростная функция фрагментации
def b(i):

	return 1.


def CalculDistribution(M,T):
	#Вычисление распределения
	
	r = random.seed()
	t = np.float64(0.0)

	increment = m.ceil( m.sqrt(M) + 2 )
	I_max =  increment

	C = np.zeros(I_max,dtype=np.int64)
	I = np.zeros(I_max,dtype=np.int64)

	NP = 1
	
	I[0] = 1
	C[0] = M

	Ac = np.zeros(I_max,dtype=np.float64)                                     # Propensity Ac[k] of coagulation of a cluster of size I[k] to size I[k]+1 
	# Склонность Ac [k] коагуляции кластера размера I [k] к размеру I [k] +1
	Tc = np.zeros(I_max,dtype=np.float64)                                     # Next time of coagulation (associate to popensity Ac[k])	
	# Следующее время коагуляции (связано с потенциалом Ac [k])
	Ab = np.zeros(I_max,dtype=np.float64)                                     # Propenstiy Ab[k] of break-up of a cluster of size I[k] to size I[k]-1
	# Склонность Ab [k] разбиения кластера размера I [k] на размер I [k] -1
	Tb = np.zeros(I_max,dtype=np.float64)                                     # Next time of break-up (associate to popensity Ab[k])			
    # следующее время разбиения
	# Init clock and propensity of each coagulation reactions.
	# инициализация часов и склонности каждой реакции коагуляции
	for i in range(0,NP):
		Ac[i] = np.float64(a(I[i])*C[i]*C[0])				
		Tc[i] = UpdateNextTime(Ac[i],Ac[i],np.float64('inf'),r,t)	

	# Init break-up of cluster size 1 to "impossible"
	# Инициализировать разбиение кластера размера 1 как "Невозможно"
	Ab[0] = np.float64('nan')
	Tb[0] = np.float64('Inf')
	
	# Init clock and propensity of each break-up reactions.
	# Инициализировать часы и склонность каждой реакции разбиения
	for i in range(1,NP):
		Ab[i] = np.float64(b(I[i])*C[i])					
		Tb[i] = UpdateNextTime(Ab[i],Ab[i],np.float64('inf'),r,t)

	# Stop the loop when a cluster of size N is created
	# Остановить цикл при создании кластера размера N
	while t < T:

		i1 = np.argmin(Tc[0:NP])
		t1 = Tc[i1]

		i2 = np.argmin(Tb[0:NP])
		t2 = Tb[i2]

		if t1<= t2:
			t  = t1
			k  = i1
			dp = 0
		
		else:
			t  = t2
			k  = i2
			dp = 1

		if t == float('inf'):
			print ('stop')
			return t

		if t <= 0:
			print ('problem negative value')
			return 

		NP = UpdateReaction(C,I,Ac,Tc,Ab,Tb,NP,k,dp,r,t)        

		if C[k] == 0 and k != 0:

			NP = ZeroReorder(C,I,Ac,Tc,Ab,Tb,NP,k)

		if NP > I_max-2:

			C  = np.concatenate((C,np.zeros(increment)))
			I  = np.concatenate((I,np.zeros(increment)))
			Ac = np.concatenate((Ac,np.zeros(increment)))
			Tc = np.concatenate((Tc,np.zeros(increment)))
			Ab = np.concatenate((Ab,np.zeros(increment)))
			Tb = np.concatenate((Tb,np.zeros(increment)))

	imax = np.max(I[0:NP])
	
	X = np.zeros(imax+1)
	Y = np.zeros(imax+1)

	for i in range(0,imax+1):	
		X[i] = i
	for k in range(0,NP):	
		Y[I[k]] = C[k]

	return X,Y


def CalculPath(M,S,T):
    #Вычисление траектории
	
	x=[]
	y=[]

	r = random.seed()
	t = np.float64(0.0)

	increment = m.ceil( m.sqrt(M) + 2 )
	I_max =  increment

	C = np.zeros(I_max,dtype=np.int64)
	I = np.zeros(I_max,dtype=np.int64)

	NP = 1
	
	I[0] = 1
	C[0] = M

	Ac = np.zeros(I_max,dtype=np.float64)                                     # Propensity Ac[k] of coagulation of a cluster of size I[k] to size I[k]+1
	# Склонность Ac [k] коагуляции кластера размера I [k] к размеру I [k] +1 
	Tc = np.zeros(I_max,dtype=np.float64)                                     # Next time of coagulation (associate to popensity Ac[k])	
	# Следующее время коагуляции (связано с потенциалом Ac [k])
	Ab = np.zeros(I_max,dtype=np.float64)                                     # Propenstiy Ab[k] of break-up of a cluster of size I[k] to size I[k]-1
    # Склонность Ab [k] разбиения кластера размера I [k] на размер I [k] -1
	Tb = np.zeros(I_max,dtype=np.float64)                                     # Next time of break-up (associate to popensity Ab[k])			
    # Следующее время разбиения
	# Init clock and propensity of each coagulation reactions.
	# инициализация часов и склонности каждой реакции коагуляции
	for i in range(0,NP):
		Ac[i] = np.float64(a(I[i])*C[i]*C[0])				
		Tc[i] = UpdateNextTime(Ac[i],Ac[i],np.float64('inf'),r,t)	

	# Init break-up of cluster size 1 to "impossible"
	# Инициализировать разбиение кластера размера 1 как "Невозможно"
	Ab[0] = np.float64('nan')
	Tb[0] = np.float64('Inf')
	
	# Init clock and propensity of each break-up reactions.
	# Инициализировать часы и склонность каждой реакции разбиения
	for i in range(1,NP):
		Ab[i] = np.float64(b(I[i])*C[i])					
		Tb[i] = UpdateNextTime(Ab[i],Ab[i],np.float64('inf'),r,t)

	# Stop the loop when the time T is reached
	# Остановить цикл при достижении времени Т
	while t < T:

		i1 = np.argmin(Tc[0:NP])
		t1 = Tc[i1]

		i2 = np.argmin(Tb[0:NP])
		t2 = Tb[i2]

		if t1<= t2:
			t  = t1
			k  = i1
			dp = 0
		
		else:
			t  = t2
			k  = i2
			dp = 1

		if t == float('inf'):
			print ('stop')
			return t

		if t <= 0:
			print ('problem negative value')
			return 

		NP = UpdateReaction(C,I,Ac,Tc,Ab,Tb,NP,k,dp,r,t)        

		if C[k] == 0 and k != 0:

			NP = ZeroReorder(C,I,Ac,Tc,Ab,Tb,NP,k)

		if NP > I_max-2:

			C  = np.concatenate((C,np.zeros(increment)))
			I  = np.concatenate((I,np.zeros(increment)))
			Ac = np.concatenate((Ac,np.zeros(increment)))
			Tc = np.concatenate((Tc,np.zeros(increment)))
			Ab = np.concatenate((Ab,np.zeros(increment)))
			Tb = np.concatenate((Tb,np.zeros(increment)))

		# output x = time, y = C_i(t)  (i=S specified size)	
		ok = 0
		for i in range(0,NP):
			if I[i] == S:		
				y = np.concatenate((y,[C[i]]))
				x = np.concatenate((x,[t]))
				ok = 1
		if ok == 0:
			y = np.concatenate((y,[0]))
			x = np.concatenate((x,[t]))

	return x,y

def CalculNucleationTime(M,N):
	#Вычисление времени зарождения
	
	r = random.seed()
	t = np.float64(0.0)

	increment = m.ceil( m.sqrt(M) + 2 )
	I_max =  increment

	C = np.zeros(I_max,dtype=np.int64)
	I = np.zeros(I_max,dtype=np.int64)

	NP = 1
	
	I[0] = 1
	C[0] = M

	Ac = np.zeros(I_max,dtype=np.float64) # Propensity Ac[k] of coagulation of a cluster of size I[k] to size I[k]+1 
    # Склонность Ac [k] коагуляции кластера размера I [k] к размеру I [k] +1 
	Tc = np.zeros(I_max,dtype=np.float64) # Next time of coagulation (associate to popensity Ac[k])
    # Следующее время коагуляции (связано с потенциалом Ac [k])
	Ab = np.zeros(I_max,dtype=np.float64) # Propenstiy Ab[k] of break-up of a cluster of size I[k] to size I[k]-1
    # Склонность Ab [k] разбиения кластера размера I [k] на размер I [k] -1
	Tb = np.zeros(I_max,dtype=np.float64) # Next time of break-up (associate to popensity Ab[k])			
    # Следующее время разбиения
	# Init clock and propensity of each coagulation reactions.
	# инициализация часов и склонности каждой реакции коагуляции
	for i in range(0,NP):
		Ac[i] = np.float64(a(I[i])*C[i]*C[0])				
		Tc[i] = UpdateNextTime(Ac[i],Ac[i],np.float64('inf'),r,t)	

	# Init break-up of cluster size 1 to "impossible"
	# Инициализировать разбиение кластера размера 1 как "Невозможно"
	Ab[0] = np.float64('nan')
	Tb[0] = np.float64('Inf')
	
	# Init clock and propensity of each break-up reactions.
	# Инициализировать часы и склонность каждой реакции разбиения
	for i in range(1,NP):
		Ab[i] = np.float64(b(I[i])*C[i])					
		Tb[i] = UpdateNextTime(Ab[i],Ab[i],np.float64('inf'),r,t)

	# Stop the loop when a cluster of size N is created
	# Остановить цикл при создании кластера размера N
	while I[NP-1] != N:

		i1 = np.argmin(Tc[0:NP])
		t1 = Tc[i1]

		i2 = np.argmin(Tb[0:NP])
		t2 = Tb[i2]

		if t1<= t2:
			t  = t1
			k  = i1
			dp = 0
		
		else:
			t  = t2
			k  = i2
			dp = 1

		if t == np.float64('inf'):
			print ('stop')
			return t

		if t <= 0:
			print ('problem negative value')
			return 

		NP = UpdateReaction(C,I,Ac,Tc,Ab,Tb,NP,k,dp,r,t)        

		if C[k] == 0 and k != 0:

			NP = ZeroReorder(C,I,Ac,Tc,Ab,Tb,NP,k)

		if NP > I_max-2:

			C  = np.concatenate((C,np.zeros(increment,dtype=np.int64)))
			I  = np.concatenate((I,np.zeros(increment,dtype=np.int64)))
			Ac = np.concatenate((Ac,np.zeros(increment,dtype=np.float64)))
			Tc = np.concatenate((Tc,np.zeros(increment,dtype=np.float64)))
			Ab = np.concatenate((Ab,np.zeros(increment,dtype=np.float64)))
			Tb = np.concatenate((Tb,np.zeros(increment,dtype=np.float64)))

	return t

def UpdateNextTime(A_new,A_old,T,r,t):

	if A_new == 0.:

		T = np.float64('inf')

	elif T == np.float64('inf'):

		U = random.random()
		T = -1./A_new*np.float64(m.log(U)) + t

	else:

		T = A_old/A_new*(T - t) + t

	return T

def UpdateCluster(C,I,NP,k,dp):

	ok = 0;
	
	# Update clusters after coagulation of index k (dp=0)
	# Обновить кластеры после коагуляции индекса k (dp = 0)
	if dp == 0:

		C[0] -= 1
		C[k] -= 1
	
		for i in range(1,NP):

			if I[i] == I[k]+1:

				C[i] += 1
				ik = i
				ok = 1

		if ok == 0:

			I[NP] = I[k] + 1
			C[NP] = 1
			ik = NP

	# Update clusters after fragmentation of index k (dp=0)
	# Обновить кластеры после фрагментации индекса k (dp = 0)
	elif dp == 1:

		C[0] += 1
		C[k] -= 1

		for i in range(0,NP):

			if I[i] == I[k]-1:

				C[i] += 1
				ik = i
				ok = 1
	        
		if ok == 0:
			I[NP] = I[k]-1
			C[NP] = 1
			ik = NP

	return ik

def UpdateReaction(C,I,Ac,Tc,Ab,Tb,NP,k,dp,r,t):

	ik = UpdateCluster(C,I,NP,k,dp)

	if ik == NP:

		NP += 1
	
		Ac[ik] = 0.
		Tc[ik] = np.float64('inf')
		Ab[ik] = 0.
		Tb[ik] = np.float64('inf')

	if dp == 0:

		Tc[k] = np.float64('inf')

		Ac_new = np.float64(a(I[0])*C[0]*(C[0]-1)/2)
		Tc[0] = UpdateNextTime(Ac_new,Ac[0],Tc[0],r,t)
		Ac[0] = Ac_new

		for i in range(1,NP):

			Ac_new = np.float64(a(I[i])*C[i]*C[0])
			Tc[i] = UpdateNextTime(Ac_new,Ac[i],Tc[i],r,t)
			Ac[i] = Ac_new

		if k != 0:

			Ab_new = np.float64(b(I[k])*C[k])
			Tb[k]  = UpdateNextTime(Ab_new,Ab[k],Tb[k],r,t)
			Ab[k]  = Ab_new

		Ab_new = np.float64(b(I[ik])*C[ik])
		Tb[ik] = UpdateNextTime(Ab_new,Ab[ik],Tb[ik],r,t)
		Ab[ik]  = Ab_new
          
	elif dp == 1:

		Ac_new = np.float64(a(I[0])*C[0]*(C[0]-1)/2)
		Tc[0] = UpdateNextTime(Ac_new,Ac[0],Tc[0],r,t)
		Ac[0] = Ac_new

		for i in range(1,NP):

			Ac_new = np.float64(a(I[i])*C[i]*C[0])
			Tc[i] = UpdateNextTime(Ac_new,Ac[i],Tc[i],r,t)
			Ac[i] = Ac_new

		if ik != 0:

			Ab_new = np.float64(b(I[ik])*C[ik])
			Tb[ik] = UpdateNextTime(Ab_new,Ab[ik],Tb[ik],r,t)
			Ab[ik] = Ab_new
		
		Tb[k] = np.float64('inf')

		Ab_new = np.float64(b(I[k])*C[k])
		Tb[k] = UpdateNextTime(Ab_new,Ab[k],Tb[k],r,t)
		Ab[k] = Ab_new

	return NP

def ZeroReorder(C,I,Ac,Tc,Ab,Tb,NP,k):

	for i in range(k,NP-1):
				
		C[i]  = C[i+1]
		I[i]  = I[i+1]
		Ac[i] = Ac[i+1]
		Tc[i] = Tc[i+1]
		Ab[i] = Ab[i+1]
		Tb[i] = Tb[i+1]					        

	return NP-1




#First exemple: Sample path of C_i for some i on [0,Tmax)
# Xp: number of paths (simulations)
# M: Number of particles

# Первый пример: Пример пути C_i для некоторого i на [0, Tmax)
# Xp: количество путей (моделирование)
# M: количество частиц
M = 1000
Xp = 10
i = 3
Tmax = 1

plt.figure(1)
plt.title('Sample path of $C_i$ with i='+str(i)+' for '+str(Xp)+'  simulations')
plt.ylabel('Normalized $C_i(t)$')
plt.xlabel('time $t$')
for n in range(0,Xp):
    [x,y] = CalculPath(M,i,Tmax)
    plt.plot(x,y,drawstyle='steps',alpha=0.75)


plt.show()

# Second exemple: Distribution of the C_i's  at time T
# Xp: number of paths (simulations)
# M: Number of particles

# Второй пример: распределение C_i во время T
# Xp: количество путей (моделирование)
# M: количество частиц
M = 1000
Xp = 10
i = 3
T = 1

plt.figure(2)
plt.title('Distirbution of $C$  at time t='+str(T)+' for '+str(Xp)+'  simulations')
plt.ylabel('Number of clusters')
plt.xlabel('size $i$')
for n in range(0,Xp):
	X,Y = CalculDistribution(M,T)
	plt.plot(X,Y,alpha=0.75)

plt.show()

#  Third exemple: Nucleation time
# T_N = inf { t > 0 : C_N(t) = 1 }
# N: size of the nucleus
# Xp: number of paths (simulations)
# M: Number of particles
# iniial condition C_1(0) = M, C_i(0)=0 for all i >1, surely.

# Третий пример: время зарождения
# T_N = inf {t> 0: C_N (t) = 1}
# N: размер ядра
# Xp: количество путей (моделирование)
# M: количество частиц
# начальное условие C_1 (0) = M, C_i (0) = 0 для всех i> 1

M = 100 #1000
Xp = 100
N = 30 #100
T = np.zeros(Xp)

for n in range(0,Xp):
    T[n] = CalculNucleationTime(M,N)

plt.figure(3)
plt.title('Nucleation time $T_N =\inf\{ t\geq 0 : C_N(t)=1\}$ for '+str(Xp)+'simulations')
plt.ylabel('Normalized number')
plt.xlabel('$T_N$')
plt.hist(T,50,facecolor='green',alpha=0.75, density = True)
plt.show()
