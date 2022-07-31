import numpy as np

N = 6

for i in range(N-1):
    for j in range(i+1,N):
        print(i,j),
    print(' ')
print(' ')


for i in range(N):
    mx = int(np.ceil(float(N-1)/2))
    if (float(N)%2 == 0.) and (i >= N/2):
        mx = int(np.floor(float(N-1)/2))

    j = i+1 - N*int(np.floor((i+1)/N + 0.5))
    cnt = 0
    while cnt < mx:
        print(i,j),
        j += 1 - N*int(np.floor((j+1)/N + 0.5))
        cnt += 1
    print(' ')

        
