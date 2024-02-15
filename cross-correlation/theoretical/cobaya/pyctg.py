from ctypes import * 
#from tqdm import tqdm

# Load the shared library
lib_correlations = CDLL("../src/libcorrelations.so")

# Define the function to calculate Ctg
ctg4py_raw=lib_correlations.ctg4py
cgg4py_raw=lib_correlations.cgg4py

# Define the function signature (argument types and return type)
ctg4py_raw.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_char_p]
ctg4py_raw.restype = c_double

cgg4py_raw.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_char_p]
cgg4py_raw.restype = c_double

lmax = 54

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
    for l in range(2, round(lmax)): #around 2-3 minutes for the whole spectrum
        #print('ctg for l=', l)
        ls.append(l)
        cl= ctg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, fname)
        ctg.append(cl)

    return ls, ctg

def cgg4py(OmegaM):

    OmegaL = 1 - OmegaM
    z0 = 0.043
    beta = 1.825
    lbda = 1.524
    h = 0.67
    bg = 1.37
    mode = 1 #mode and ncalls don't change the results, bur are still needed
    ncalls = 200000
    pkfname = "../tables/pk_3dmatter.dat"
    fname   = c_char_p(pkfname.encode("ascii"))
    ls=[]
    cgg = []
    for l in range(2, round(lmax)): #about 12 seconds per point, ~6mins for 54 points
        #print("cgg for l=", l)
        ls.append(l)
        cl= cgg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, fname)
        cgg.append(cl)

    return ls, cgg

###Testing

if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    ls,ctg=ctg4py(0.3)
    print('ls=', ls)
    print('ctg(0.3)=',ctg)

    ls,cgg=cgg4py(0.3)
    print('ls=', ls)
    print('cgg(0.3)=',cgg)


    plt.figure()
    plt.plot(ls,ctg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{tg}$')
    plt.xscale('log')
    plt.savefig("pyctg_full_test.png")

    plt.figure()
    plt.plot(ls, cgg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{gg}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("pycgg_full_test.png")

