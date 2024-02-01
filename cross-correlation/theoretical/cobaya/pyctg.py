from ctypes import * 

# Load the shared library
lib_correlations = CDLL("../src/libcorrelations.so")

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
    ncalls = 200000
    pkfname = "../tables/pk_3dmatter.dat"
    fname   = c_char_p(pkfname.encode("ascii"))
    ls=[]
    ctg = []
    for l in range(2, round(lmax)):
        ls.append(l)
        cl= ctg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, fname)
        ctg.append(cl)

    return ls, ctg
    

#Still need to implement Cgg?

###Testing

if __name__=='__main__':
    import matplotlib.pyplot as plt

    ls,ctg=ctg4py(0.3)
    print('ls=', ls)
    print('ctg(0.3)=',ctg)

    plt.figure()
    plt.plot(ls,ctg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_\ell^{tg}$')
    plt.xscale('log')
    plt.savefig("pyctg_full_test.png")

