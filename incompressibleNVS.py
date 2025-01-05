import numpy as np
import matplotlib.pyplot as plt

X_CELL_NUM = 20
Y_CELL_NUM = 20

DIST = 1
DELTA_T = 0.02

MU = 0.5
RO = 20

#np.seterr(invalid='ignore')

# ------------------------------------------------------------------------laplacien as matrix operation------------------------------------------------------------------------
def map_2d_to_1d(x, y, size_y):
    return y+x*size_y

def calc_stencil_matrix_laplacian(size_x, size_y):
    # this almost broke my brain
    # you can generalize this way better for arbitrary stencils i think but for now im happy
    # this maps the common laplacien stencil 010,141,010 to a matrix A to use as a linear system of eq A*b = c for implicitly solving the poisson equation for pressure
    finite_difference_matrix = np.diag(size_x*size_y*[0])
    
    for x_c in range(0, size_x):
        for y_c in range(0, size_y):
            finite_difference_matrix[map_2d_to_1d(x_c, y_c, size_y)][map_2d_to_1d(x_c, y_c, size_y)] = -4
            if y_c > 0:
                finite_difference_matrix[map_2d_to_1d(x_c, y_c-1, size_y)][map_2d_to_1d(x_c, y_c, size_y)] = 1
            if y_c < size_y-1:
                finite_difference_matrix[map_2d_to_1d(x_c, y_c+1, size_y)][map_2d_to_1d(x_c, y_c, size_y)] = 1
            if x_c > 0:
                finite_difference_matrix[map_2d_to_1d(x_c-1, y_c, size_y)][map_2d_to_1d(x_c, y_c, size_y)] = 1
            if x_c < size_x-1:
                finite_difference_matrix[map_2d_to_1d(x_c+1, y_c, size_y)][map_2d_to_1d(x_c, y_c, size_y)] = 1
                 
    return finite_difference_matrix 
    
LAPLACIAN_MATRIX = calc_stencil_matrix_laplacian(X_CELL_NUM, Y_CELL_NUM)
print(calc_stencil_matrix_laplacian(3, 3))

def mat_to_vec(m):
    return m.flatten()
    
def vec_to_mat(v, size_x, size_y):
    return v.reshape(size_x, size_y)    

#------------------------------------------------------------------------finite difference helper functions--------------------------------------------------------------------
# I'm aware that this code can be much faster with numpy matrices, this is just for experimentation 

def central_fin_dif_x(scalar_field, size_x, size_y, h):
    result = scalar_field.copy()
    for x in range(1, size_x-1):
        for y in range(1, size_y-1):
            result[x][y] = (scalar_field[x+1][y] - scalar_field[x-1][y])/(2*h)
    return result
    
def central_fin_dif_y(scalar_field, size_x, size_y, h):
    result = scalar_field.copy()
    for x in range(1, size_x-1):
        for y in range(1, size_y-1):
            result[x][y] = (scalar_field[x][y+1] - scalar_field[x][y-1])/(2*h)
    return result
    
def divergence(scalar_field_x, scalar_field_y, size_x, size_y, h):
    a = central_fin_dif_x(scalar_field_x, size_x, size_y, h)
    b = central_fin_dif_y(scalar_field_y, size_x, size_y, h)
    return a + b

# doesnt work apparently
def laplacian(scalar_field, size_x, size_y, h):
    a = central_fin_dif_x(central_fin_dif_x(scalar_field, size_x, size_y, h), size_x, size_y, h)
    b = central_fin_dif_y(central_fin_dif_y(scalar_field, size_x, size_y, h), size_x, size_y, h)
    return a + b
    
def laplacian2(scalar_field, size_x, size_y, h):
    result_a = scalar_field.copy()
    result_b = scalar_field.copy()
    for x in range(1, size_x-1):
        for y in range(1, size_y-1):
            result_a[x][y] = (scalar_field[x+1][y] - 2*scalar_field[x][y] + scalar_field[x-1][y])/(h*h)
    for x in range(1, size_x-1):
        for y in range(1, size_y-1):
            result_b[x][y] = (scalar_field[x][y+1] - 2*scalar_field[x][y] + scalar_field[x][y-1])/(h*h)
    return result_a + result_b
    
def convective_derivative(u,v,prop, size_x, size_y, h):
    return u*central_fin_dif_x(prop, size_x, size_y, h) + v*central_fin_dif_y(prop, size_x, size_y, h)
    
def solve_poisson(div_H, size_x, size_y):      
    p = vec_to_mat(np.linalg.solve(LAPLACIAN_MATRIX, mat_to_vec(div_H)), size_x, size_y)
    return p
    
def solve_poisson2(div_H, u, v, ro, dt):
    pass
    
#------------------------------------------------------------------------------main loop and plotting------------------------------------------------------------

def main():
    # u, v should be divergence free
    u = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    v = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    p = np.zeros((X_CELL_NUM, Y_CELL_NUM))
    
    # plt init
    x, y = np.meshgrid(np.linspace(0, X_CELL_NUM, X_CELL_NUM),np.linspace(0, Y_CELL_NUM, Y_CELL_NUM)) 
    
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4)
    
    # for i in range(5):
    #    v[5][10] = 100
    #u[10][30] = -100
    #p[9][31] = 50
    
    pressure_history = []

    
    for t in range(100000):
        # calculate H_i
        H_u = MU*laplacian2(u, X_CELL_NUM, Y_CELL_NUM, DIST)-convective_derivative(u, v, u, X_CELL_NUM, Y_CELL_NUM, DIST)
        H_v = MU*laplacian2(v, X_CELL_NUM, Y_CELL_NUM, DIST)-convective_derivative(u, v, v, X_CELL_NUM, Y_CELL_NUM, DIST)
        
        div_H = divergence(H_u, H_v, X_CELL_NUM, Y_CELL_NUM, DIST)
        
        # solve poisson equation to obtain pressure field
        p = solve_poisson(div_H, X_CELL_NUM, Y_CELL_NUM)
        
        max_pressure_grad = p.max()-p.min()
        pressure_history.append(max_pressure_grad)
        print("-------------------------------")
        print(f"Iteration: {t}")
        print(f"max pressure gradient: {max_pressure_grad}")
        print(f"mean pressure: {p.mean()}")
        print(f"mean u: {u.mean()}")
        print(f"mean v: {v.mean()}")
        print(f"max divergence of velocity: {divergence(u,v, X_CELL_NUM, Y_CELL_NUM, DIST).max()}")
        
        u+=(-central_fin_dif_x(p, X_CELL_NUM, Y_CELL_NUM, DIST)*(1/RO) + H_u)*DELTA_T
        v+=(-central_fin_dif_y(p, X_CELL_NUM, Y_CELL_NUM, DIST)*(1/RO) + H_v)*DELTA_T
        
        for x in range(X_CELL_NUM):
            u[x][0] = 0
            u[x][-1] = 0
            v[x][0] = 0
            v[x][-1] = 0
        for y in range(Y_CELL_NUM):
            u[0][y] = 0
            u[-1][y] = 0
            v[0][y] = 10
            v[-1][y] = 0
        
        ax1.clear()
        ax1.set_title("Pressure")
        ax1.imshow(p)
        ax2.clear()
        ax2.set_title("u")
        ax2.imshow(u)
        ax3.clear()
        ax3.set_title("v")
        ax3.imshow(v)
        ax4.clear()
        ax4.set_title("max pressure gradient")
        ax4.plot(pressure_history)
        plt.pause(0.01)
        
if __name__ == "__main__":
    main()
    
    


