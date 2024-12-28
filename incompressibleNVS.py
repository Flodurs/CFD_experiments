import numpy as np
import matplotlib.pyplot as plt

X_CELL_NUM = 100
Y_CELL_NUM = 100

DIST = 0.01
DELTA_T = 0.05

MU = 1
RO = 1



# I'm aware that this code can be much faster with numpy matrices, this is just for experimentation 

def central_fin_dif_x(scalar_field, size_x, size_y, h):
    result = scalar_field.copy()
    for x in range(1, size_x-1):
        for y in range(0, size_y-1):
            result[x][y] = (scalar_field[x+1][y] - scalar_field[x-1][y])/(2*h)
    return result
    
def central_fin_dif_y(scalar_field, size_x, size_y, h):
    result = scalar_field.copy()
    for x in range(0, size_x-1):
        for y in range(0, size_y-1):
            result[x][y] = (scalar_field[x][y+1] - scalar_field[x][y-1])/(2*h)
    return result
        
def test_central_fin_dif():
    A = np.diag([1,2,3,4,5])
    print(A)
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print(central_fin_dif_x(A,5,5,1))
    print("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY")
    print(central_fin_dif_y(A,5,5,1))
    
def divergence(scalar_field_x, scalar_field_y, size_x, size_y, h):
    a = central_fin_dif_x(scalar_field_x, size_x, size_y, h)
    b = central_fin_dif_y(scalar_field_y, size_x, size_y, h)
    return a + b
    
def laplacian(scalar_field_x, scalar_field_y, size_x, size_y, h):
    a = central_fin_dif_x(central_fin_dif_x(scalar_field_x, size_x, size_y, h), size_x, size_y, h)
    b = central_fin_dif_y(central_fin_dif_y(scalar_field_y, size_x, size_y, h), size_x, size_y, h)
    return a + b
    
def solve_poisson(dif_H):
    pass
    


def main():
    # u, v should be divergence free
    u = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    v = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    p = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    
    # plt init
    x, y = np.meshgrid(np.linspace(0, X_CELL_NUM, X_CELL_NUM),np.linspace(0, Y_CELL_NUM, Y_CELL_NUM)) 
    
    fig, ax = plt.subplots()
    
    u[20][30] = 1
    
    for t in range(100):
        # calculate H_i
        H_u = MU*laplacian(u, v, X_CELL_NUM, Y_CELL_NUM, DIST)*1/RO-np.multiply(divergence(u, v, X_CELL_NUM, Y_CELL_NUM, DIST), u)
        H_v = MU*laplacian(u, v, X_CELL_NUM, Y_CELL_NUM, DIST)*1/RO-np.multiply(divergence(u, v, X_CELL_NUM, Y_CELL_NUM, DIST), v)
        div_H = divergence(H_u, H_v, X_CELL_NUM, Y_CELL_NUM, DIST)
        u+=H_u*DELTA_T
        v+=H_v*DELTA_T
        ax.clear()
        ax.imshow(u)
        plt.pause(0.01)
        
        
main()
    
    


