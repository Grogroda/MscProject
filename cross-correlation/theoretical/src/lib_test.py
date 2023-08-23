from ctypes import * 

# Load the shared library
libc = CDLL("./libcorrelations.so")

# Define the function signature (argument types and return type)
ctg=libc._Z3ctgddidddddii
print("type=", type(ctg))
ctg.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int]
ctg.restype = c_double

OmegaL = 0.7
Omegam = 0.3
l = 11
z0 = 0.15
beta = 3.09 
lbda = 4.94
h = 0.67
bg = 1.0
mode = 1
ncalls = 1000
print("Starting calculation")
#result = libc._Z3ctgddidddddii(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls)
result = ctg(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls)
print("Result=", result)
