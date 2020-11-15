import numpy as np

lx = 5000.0
nx = [40, 50, 60, 80, 100]
dt1=5*60
nx1=50
dx1 = lx/nx1

dt = np.zeros(6)
btr_dt = np.zeros(6)
dx = np.zeros(6)

print('j    nx   dx   nx*dx   dt   btr_dt')
for j in range(len(nx)):
    dx[j] = lx / nx[j]
    dt[j] = dt1 *dx[j] / dx1
    btr_dt[j] = dt[j] / 20.0
    print(j,nx[j],dx[j],nx[j]*dx[j],dt[j], btr_dt[j])


