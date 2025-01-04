import numpy as np


def map_2d_to_1d(x,y,size_y):
    return y+x*size_y

def calc_stencil_matrix_laplacian(size_x, size_y):
    # this almost broke my brain
    # you can generalize this way better for arbitrary stencils i think but for now im happy
    finite_difference_matrix = np.diag(size_x*size_y*[0])
    for x_c in range(size_x):
        for y_c in range(size_y):
            print(map_2d_to_1d(x_c, y_c, size_y))
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

def create_laplacien_matrix(size_x, size_y):
    m = np.diag([-4]*size_x*size_y) 
    l = np.diag([1]*(size_x*size_y-size_y), k = size_y)
    u = np.diag([1]*(size_x*size_y-size_x), k = -size_x)
    
    loechrig_up = np.zeros(size_x*size_y-1)
    for i in range(len(loechrig_up)):
        if (i+5)%5 == 0 :
            loechrig_up[i:i+5] = 1
     
    loechrig_low = np.zeros(size_x*size_y-1)
    for i in range(len(loechrig_low)):
        if (i+5)%5 == 0 :
            loechrig_low[i:i+6] = 1  
            
            
    ll = np.diag(loechrig_low, k=-1)
    lu = np.diag(loechrig_up, k=1)
    return m+u
    
def mat_to_vec(m):
    return m.flatten()
    
def vec_to_mat(v, size_x, size_y):
    return v.reshape(size_x, size_y)
    
    
    
m = calc_stencil_matrix_laplacian(5,6)
print(m)
d = np.zeros((5,6))

d[0][0] = 3
d[4][5] = 3



A = np.matmul(m,mat_to_vec(d))

A = vec_to_mat(A,5,6)

print(A)