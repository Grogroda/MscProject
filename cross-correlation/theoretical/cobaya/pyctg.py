from ctypes import * 

# Load the shared library
lib_correlations = CDLL("./libcorrelations.so")

# Define the function to calculate Ctg
ctg4py_raw=lib_correlations.ctg4py

# Define the function signature (argument types and return type)
ctg4py_raw.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_char_p]
ctg4py_raw.restype = c_double

lmax = 6

def ctg4py(OmegaM):
    OmegaL = 1 - OmegaM
    z0 = 0.043
    beta = 1.825 
    lbda = 1.524
    h = 0.67
    bg = 1.37
    mode = 1
    ncalls = 50000
    pkfname = "../tables/pk_3dmatter.dat"
    fname   = c_char_p(pkfname.encode("ascii"))
    ctg = []
    for l in range(2, round(lmax)):
        cl= ctg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, fname)
        ctg.append(cl)

    return ctg
    

#Still need to implement Cgg?

###Testing

if __name__=='__main__':
    import matplotlib.pyplot as plt

    ls=[2+i for i in range(lmax-2)]
    ctg=ctg4py(0.3)
    print('ctg(0.3)=',ctg)

    plt.figure()
    plt.plot(ls,ctg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_\ell^{tg}$')
    plt.xscale('log')
    plt.savefig("pyctg_full_test.png")

"""
if __name__=='__main__':
    #Testing area:

    # Define test parameters
    OmegaL = 0.7
    Omegam = 0.3
    l = 9 
    z0 = 0.15
    beta = 3.09 
    lbda = 4.94
    h = 0.67
    bg = 1.0
    mode = 1
    ncalls = 1000
    pkfname = "../tables/pk_3dmatter.dat"
    fname=c_char_p(pkfname.encode("ascii"))

    # Test shared library:
    print("Starting calculation")
    for l in [2,8,12]:
        result = ctg4py(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls, fname)
        print("Result=", result)
"""
