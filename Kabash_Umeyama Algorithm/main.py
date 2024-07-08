import numpy as np
import math

def read_points(filename):
    f = open(filename, 'r')
    lines = f.readlines()

    points = np.array([[float(val) for val in line.strip().split()] for line in lines])
    f.close()
    #print(points)
    return points

def read_correspondences(filename):
    correspondences = []
    f = open(filename, 'r')
    for line in f:
        point1, point2 = map(float, line.strip().split())
        correspondences.append((point1, point2))
    f.close()
    #print(correspondences)
    return correspondences

def determinant(A):
    boyut = np.shape(A)[0]
    if boyut == 2:
       return A[0][0]*A[1][1] - A[0][1]*A[1][0]
    else:
        sign = 1
        cumulative_sum = 0
        col_num = 0
        while col_num < boyut:
            sub_matrix = np.zeros((boyut-1,boyut-1))
            for i in range(1, boyut):
                flag = False
                for j in range(boyut):
                    if(j == col_num):
                        flag = True
                    else:
                        if flag:
                            sub_matrix[i-1][j-1] = A[i][j]
                        else:
                            sub_matrix[i-1][j] = A[i][j]
            det = determinant(sub_matrix)
            cumulative_sum += A[0,col_num] * det * sign
            sign = -sign
            col_num += 1   
        return cumulative_sum

#A = np.array([[1,2,6],[2,5,7],[1,2,7]])#,[1,1,1],[1,0,9]
#print(determinant(A))

def norm(e):
    # e: 1x3 vector
    return math.sqrt(e[0]**2 + e[1]**2 + e[2]**2)

def qr(A):
    #A nx3 matrix
    At = A.T # At: 3xn
    e1 = At[0] #1xn
    e2 = At[1] - ((At[1].T)@e1)/((e1.T)@e1)*e1
    e3 = At[2] - ((At[2].T)@e1)/((e1.T)@e1)*e1 - ((At[2].T)@e2)/((e2.T)@e2)*e2
    
    #print(e1)
    #print(np.linalg.norm(e1))
    #print(norm(e1))
    
    e1 = e1/norm(e1)
    #print(e1)
    e2 = e2/norm(e2)
    e3 = e3/norm(e3)

    Q = np.array([e1, e2, e3]).T
    R = (Q.T)@A
    return Q, R
    
#print(np.linalg.qr(A))
#print(qr(A))
#qr(A)

def eig(A):

    A_copy = np.copy(A)
    eigenvectors = np.eye(np.shape(A)[0])
    
    for i in range(200):
        
        #Q, R = np.linalg.qr(A_copy)
        Q, R = qr(A_copy)
        A_copy = R@Q
        
        off_diagonal = np.abs(A_copy - np.diag(np.diag(A_copy))).sum()
        if off_diagonal < 1e-6:
            break
        eigenvectors = np.dot(eigenvectors, Q)
        
    eigenvalues = np.diag(A_copy)
    
    return eigenvalues, eigenvectors

#print(eig(A))
    
def svd(A):
    ATA = (A.T)@A
    AAT = A@(A.T)
    eigenvalues_ATA, U = eig(ATA)
    eigenvalues_AAT, V = eig(AAT)

    S = np.sqrt(eigenvalues_ATA)
    idx = np.argsort(-S)
    S = S[idx]
    U = U[:, idx]
    V = V[:, idx]

    Vt = V.T

    return U, S, Vt

#print(svd(A))
#A2 = np.zeros((np.shape(A)[0],3))
#print(A2)

def kabash_umeyama(Q, P, correspondences): 
    Q_matched = np.zeros((np.shape(correspondences)[0],3))
    count = 0
    for point1, i in correspondences:
        Q_matched[count] = Q[int(point1)]
        count += 1
        
    #print(Q_matched)
    
    P_matched = np.zeros((np.shape(correspondences)[0],3))
    count = 0
    for i, point2 in correspondences:
        P_matched[count] = P[int(point2)]
        count += 1
    
    #print(P_matched)
    
    Q_centered = Q_matched - np.mean(Q_matched, axis=0)
    P_centered = P_matched - np.mean(P_matched, axis=0)

    H = (Q_centered.T)@(P_centered)

    U, S, Vt = svd(H)
    R = (Vt.T)@(U.T)
    if determinant(R) < 0:
        Vt[-1] *= -1
        R = (Vt.T)@(U.T)
        
    t = np.mean(P_matched, axis=0) - R@(np.mean(Q_matched, axis=0))

    return R, t

mat1 = read_points("mat1.txt")
mat2 = read_points("mat2.txt")
correspondences = read_correspondences("correspondences.txt")

R, t = kabash_umeyama(mat1, mat2, correspondences)

mat2_transformed = mat2@(R.T) + t

merged_points = np.concatenate((mat1, mat2_transformed))

np.savetxt("rotation_mat.txt", R.T, fmt="%0.6f")
np.savetxt("translation_vec.txt", t, fmt="%0.6f")
np.savetxt("merged.txt", merged_points, fmt="%0.6f")
