import numpy as np
import matplotlib.pyplot as plt
import incompressibleNVS

def test_numpy_array():
    A = np.diag((1,2,3,4,5))
    A[1][3] = 3
    print(A)
    plt.imshow(A)
    plt.show()
    
def test_laplacian():
    A = np.zeros((5,5))     
    A[2][3] = 1
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.imshow(A)
    ax2.imshow(incompressibleNVS.laplacian2(A,5,5,1))
    plt.show()
        
test_laplacian()