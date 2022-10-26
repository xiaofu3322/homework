import numpy as np
import matplotlib.pyplot as plt





def Lagrange(xk, yk, x):
        n = len(xk)
        lk = []    #后面的项
        y = 0
        for k in range(n):
            zi = 1
            mu = 1
            for j in range(n):
                if j != k:
                    zi *= (x - xk[j])       #分子
                    mu *= (xk[k] - xk[j])   #分母
            lk.append(zi / mu)

        for i in range(n):
            y += lk[i] * yk[i]
        return y

def newton(xk,yk,x):
    n = len(xk)
    y = yk[0]                       #先把y0加入最后的项中
    A = np.zeros([n,n])
    for i in range(0,n):
        A[i,0]=yk[i]

    sum1 = 1   #x的多项式
    for i in range(1,n):
        sum1 = sum1*(x-xk[i-1])
        for j in range(i, n):
            A[j,i] = (A[j,i-1]-A[j-1,i-1])/(xk[j]-xk[j-i])
        y+=sum1*A[1,1]  #对角
    return y

xk = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
yk = [1/26, 1/17, 0.1, 0.2, 0.5, 1, 0.5, 0.2, 0.1, 1/17, 1/26]

xs=np.linspace(np.min(xk),np.max(xk),1000,endpoint=True)
yl=[]
yn=[]
for x in xs:
    yl.append(Lagrange(xk,yk,x))
    yn.append(newton(xk,yk,x))

lines1 = plt.plot(xs,yl)
lines2 = plt.plot(xs,yn)

plt.setp(lines1, color='r', linewidth=4.0)
plt.setp(lines2, color='y', linewidth=2.0)

plt.show()

print(lines)
