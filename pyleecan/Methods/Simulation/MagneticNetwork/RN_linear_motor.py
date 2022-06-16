# %%
# %%
# Machine's characteristics





import time
from threadpoolctl import threadpool_info
import os

threads = "1"
os.environ["OMP_NUM_THREADS"] = threads  # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = threads  # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = threads  # export MKL_NUM_THREADS=6
# export VECLIB_MAXIMUM_THREADS=4
os.environ["VECLIB_MAXIMUM_THREADS"] = threads
os.environ["NUMEXPR_NUM_THREADS"] = threads  # export NUMEXPR_NUM_THREADS=6
# export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=2"
os.environ["MKL_DOMAIN_NUM_THREADS"] = "MKL_DOMAIN_ALL=1"


os.environ["MKL_DYNAMIC"] = "FALSE"
os.environ["OMP_DYNAMIC"] = "FALSE"

import numpy as np

print("First call:", threadpool_info())
print(np.__config__.show())

#Choose beetween chol-mod and scipy sparse linear solver
Have_cholmod=False
try:
    from sksparse.cholmod import cholesky,analyze
    Have_cholmod=True
except:
    from scipy.sparse.linalg import spsolve, bicg, cg

from scipy.sparse import coo_matrix
from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import dia_matrix

import matplotlib.pyplot as plt



# %%
# Machine's characteristics
tp = 60e-3        # pole pitch (m)
tm = 55e-3        # PM length in x direction (m)
hm = 10e-3        # PM height in y direction (m)
e = 1e-3         # Air-gap thickness (m)
hst = 30e-3        # Stator total height (m)
hs = 20e-3        # Slot height (m)
hmbi = 10e-3        # Moving armature height (moving back iron height)
ws = 10e-3        # Slot opening (m)
ts = 2*ws         # Slot pitch (m)
la = 1            # Active length (m)
Br = 1.2          # PM remanent induction (residual induction) (T)
mu0 = np.pi*4e-7    # Permeability of vacuum (H/m)
mur1 = 1            # Relative permeability of air
mur2 = 7500         # Relative permeability of statot iron
mur3 = 7500

# %%


size_x =61
size_y = 52

x_min = 0
x_max = tp


y_min = 0
y_max = e+hst+hm+hmbi


h_x = (x_max-x_min)/(size_x-1)
h_y = (y_max-y_min)/(size_y-1)


# %%
# % Number of elements in the stator armature
m0s = round((ts-ws)/2/h_x)  # Number of elements in half a tooth in x direction
m1s = round(ws/h_x)
# Number of elements in the stator back iron in y direction
p0s = round((hst-hs)/h_y)
m = 12*m0s               # Total number of elements of the stator in x direction
p = 3*p0s                # Total number of element of stator in y direction

# Number of elements in the moving armature (the air-gap is supposed to be part of the moving armature)
# Number of elements in half the air-gap between two adjacent PM in x direction
m0m = round((tp-tm)/2/h_x)
m1m = round(tm/2/h_x)
p0m = round(e/h_y)          # Number of elements in the air-gap in y direction
# Number of elements in the moving armature iron in y direction
p0 = round(hmbi/h_y)
# Number of elements in the magnetic air-gap (hm + e) in y direction
p1 = round((hm+e)/h_y)


# %%


x = np.linspace(x_min, x_max, size_x)
y = np.linspace(y_min, y_max, size_y)


# x=(x_max-x_min)*np.sort(np.random.random_sample(size_x))+x_min
# y=(y_max-y_min)*np.sort(np.random.random_sample(size_y))+y_min


x_dual = (x[1:]+x[:size_x-1])/2
y_dual = (y[1:]+y[:size_y-1])/2





list_materials = ["bob1", "bob2", "bob3", "air", "iron1", "PM", "iron3"]

permeabilty_materials = np.array([1, 1, 1, 1, 7500, 1, 7500])




def geometry(size_x, size_y,pos_pm):
    m = size_y-1
    n = size_x-1
    nn = m*n

    cells_materials = np.zeros(nn, dtype=np.uint16)

    mask_magnet=np.zeros(nn,dtype=np.bool_)
    mask_magnet[n*p1-n:n*(p1+p0-p0m)+n]=True
    ### Geometry assembly
    for i in range(m):
        if m-p0s <= i:
            for j in range(n):
                num_elem = n*i+j
                cells_materials[num_elem] = 5

        elif (p0+p1 <= i < n-p0s):
            for j in range(n):
                num_elem = n*i+j
                if (m0s <= j < m0s+m1s):
                    cells_materials[num_elem] = 1
                elif (3*m0s+m1s <= j < 3*m0s+2*m1s):
                    cells_materials[num_elem] = 2
                elif (5*m0s+2*m1s <= j < 5*m0s+3*m1s):
                    cells_materials[num_elem] = 3
                else:
                    cells_materials[num_elem] = 5
        
        elif p0+p1-p0m <= i < p0+p1:
            for j in range(n):
                num_elem = n*i+j
                cells_materials[num_elem] = 4
        ##
        elif p1 <= i < p1+p0-p0m:
            for j in range(n):
                num_elem = n*i+j
                if pos_pm+2*m1m>=n:
                    if (pos_pm+2*m1m)%n<j<=pos_pm%n:
                        cells_materials[num_elem] = 4
                    else:
                        cells_materials[num_elem] = 6
                else:
                    if pos_pm <= j < (pos_pm+2*m1m):
                        cells_materials[num_elem] = 6
                    else:
                        cells_materials[num_elem] = 4
        ##
        elif i < p1:
            for j in range(n):
                num_elem = n*i+j
                cells_materials[num_elem] = 7
        else:
            print("Wrong geometry")
    

    return cells_materials








def linear_model(size_x, size_y, x, y, x_dual, y_dual, pos,BC):
    t0 = time.perf_counter()
    # Initialyze nature of elements
    list_geometry = geometry(size_x, size_y, pos)

    # Initialyze mesh
    #Num_Unknowns -> numerotation in the matrix
    #list_elem -> the list of the elements (cells)
    #BC_list -> Boudary condition
    
    list_coord=init_point(size_x,size_y,x,y)
    
    permeability_cell=init_permeabilty_cell(size_x,size_y,permeabilty_materials,list_geometry)
    
    list_elem=init_cell(size_x, size_y,x,y)
    
    BC_list,Periodic_point = init_mesh_BC(size_x, size_y,x,y,BC)
    
    Num_Unknowns=numeroting_unknows(list_elem,BC_list,Periodic_point)
    
    
    t1 = time.perf_counter()
    print("Assembly geometry:", np.round(t1-t0, 5), "secondes")
    save_mesh(list_geometry,Num_Unknowns,list_elem,x,y,BC_list)
    t2 = time.perf_counter()
    print("Save mesh:", np.round(t2-t1, 5), "secondes")

    # Assembly all matrice
    reluc_list=init_reluc(list_elem,list_coord)
    M_csr = assembly(reluc_list, Num_Unknowns, list_elem, permeability_cell,BC_list)
    
    t3 = time.perf_counter()
    print("Assembly matrix", np.round(t3-t2, 5), "secondes")

    # Assembly RHS
    E = right_member_assembly(list_geometry, Num_Unknowns, list_elem)
    t4 = time.perf_counter()
    
    print("Assembly vector:", np.round(t4-t3, 5), "secondes")
    print("Total :", np.round(t4-t2, 5))
    
    #Compute Solution
    t3= time.perf_counter()
    if Have_cholmod:
        
        # Compute the solution
        factor = cholesky(M_csr.tocsc())
        F = factor(E)
        t4 = time.perf_counter()
        print("Time to solve (CholMod):", np.round(t4-t3, 5),
              "secondes, res:", np.linalg.norm(M_csr@F-E, ord=2))
    else:
        F = spsolve(M_csr, E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)
        t4 = time.perf_counter()
        print("Time to solve (direct,UMF):", np.round(t4-t3, 5),
              "secondes, res:", np.linalg.norm(M_csr@F-E, ord=2))
    
    F=add_BC_to_F(F,size_x,size_y,Num_Unknowns, list_elem,BC_list)
    
    return F, list_geometry, Num_Unknowns, list_elem, permeability_cell,list_coord


def mu_non_linear(Phi, permeability_cell, list_geometry,Num_Unknowns,list_elem,BC_list):
    # Compute staured mu
    B_sat = 1.8
    materials = 2
    mu_r0 = 7500
    Phi=add_BC_to_F(Phi,size_x,size_y, Num_Unknowns,list_elem,BC_list)
    Bx, By = compute_B_square(Phi, size_x, size_y, list_geometry)
    B_norm = np.sqrt(Bx*Bx+By*By).flatten()

    H = B_norm/(mu0*permeability_cell)
    mask = list_geometry == materials
    print("H=", np.round(H[mask].max(),5), "A/m, mu=", np.round(permeability_cell[mask].min(),5),"H/m")

    permeability_cell[mask] = 1+(2*B_sat/(np.pi*mu0*H[mask])) * \
        np.arctan(np.pi*(mu_r0-1)*mu0*H[mask]/(2*B_sat))
    # print(permeability_cell[mask])
    return permeability_cell


def non_linear_model(size_x, size_y, x, y, x_dual, y_dual, pos, sol0,BC):
    t0 = time.perf_counter()
    # Initialyze nature of elements
    list_geometry = geometry(size_x, size_y, x, y, x_dual, y_dual, pos)

    # Initialyze mesh
    list_coord=init_point(size_x,size_y,x,y)
    
    permeability_cell=init_permeabilty_cell(size_x,size_y,permeabilty_materials,list_geometry)
    
    list_elem=init_cell(size_x, size_y,x,y)
    
    BC_list,Periodic_point = init_mesh_BC(size_x, size_y,x,y,BC)
    
    Num_Unknowns=numeroting_unknows(list_elem,BC_list,Periodic_point)
    
    t1 = time.perf_counter()
    print("Assembly geometry:", np.round(t1-t0, 5), "secondes")
    # Assembly each area of the matrix
    t2 = time.perf_counter()
    Matrices = []
    print(permeabilty_materials)
    for i in range(1, permeabilty_materials.size+1):
        M_csr = assembly_one_area(i, list_materials, permeabilty_materials, size_x, size_y,
                                  x, y, x_dual, y_dual, list_geometry, Num_Unknowns, list_elem, permeability_cell,BC_list)
        Matrices.append(M_csr)

    t3 = time.perf_counter()
    print("Assembly matrix", np.round(t3-t2, 5), "secondes")

    # Assembly RHS
    E = right_member_assembly(list_materials, permeabilty_materials, size_x, size_y,
                              x, y, x_dual, y_dual, list_geometry, Num_Unknowns, list_elem, permeability_cell)
    t4 = time.perf_counter()
    print("Assembly vector:", np.round(t4-t3, 5), "secondes")
    print("Total :", np.round(t4-t2, 5))
    t12 = time.perf_counter()
    M_csr= Matrices[0]
    # Togather area of the matrix
    for i in range(permeabilty_materials.size-1):
        M_csr += M_csr + Matrices[i+1]
    M2 = Matrices[0]+Matrices[2]

    # Initialize non linear loop
    
    
    if Have_cholmod:
        M_csc=M_csr.tocsc()
        #Symbolic Cholesky factorisation
        factor=analyze(M_csc,mode="simplicial",ordering_method="nesdis")
        #Cholesky factorisation
        factor.cholesky_inplace(M_csc)
        F = factor(E)
    else:
        F = spsolve(M_csr, E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)

    
    t1 = time.perf_counter()
    old_permaebility = permeability_cell.copy()
    phi_old = F.copy()
    i = 0
    eps = 10**-4
    criteria = eps+1
    # Non linear loop
    
    while criteria > eps and i < 500:
        print("loop n°", i+1)
        # Update permeability
        permeability_cell = mu_non_linear(F, permeability_cell, list_geometry,Num_Unknowns,list_elem,BC_list)
        # Update matrix
        M3 = assembly_one_area(2, list_materials, permeabilty_materials, size_x, size_y,
                               x, y, x_dual, y_dual, list_geometry, Num_Unknowns, list_elem, permeability_cell,BC_list)
        if Have_cholmod:
           # Compute cholesky decomposition, solve the system, compute solution
           factor.cholesky_inplace((M2+M3).tocsc())
           F = factor(E)
        else:
           F = spsolve((M2+M3), E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True) 
        
       
        # Stoping criterion
        criteria = np.linalg.norm(F-phi_old, ord=2)/np.linalg.norm(F, ord=2)
        criteria2= np.linalg.norm(old_permaebility-permeability_cell, ord=2)/np.linalg.norm(permeability_cell, ord=2)
        # Updtae old solution
        phi_old = F.copy()
        old_permaebility = permeability_cell.copy()
        i += 1
        print("eps flux:", np.round(criteria,7),"eps permeability", np.round(criteria2,9) ,"\n")
    
    t2 = time.perf_counter()
    print("Time to solve non-linear:", np.round(t2-t1, 5),
          "secondes, res:", np.linalg.norm((M2+M3)@F-E, ord=2))
    
    F=add_BC_to_F(F,size_x,size_y,Num_Unknowns, list_elem,BC_list)
    return F, list_geometry




