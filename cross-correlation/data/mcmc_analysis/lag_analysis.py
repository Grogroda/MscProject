import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ac(l, sample):
    '''
    l is an int for the value of the lag and sample is the list of sampled values for the variable of interest. In my case, sample will be either samples[Ctt], samples[Cgg] or samples[Ctg] for each multipole number.
    '''
    num=0
    den=0

    N=len(sample)
    
    barx=1/N*sum(sample)

    for i in range(N): 
            den+=(sample[i]-barx)**2
            if i>=l:
                num+=(sample[i]-barx)*(sample[i-l]-barx)

    num=num*1/(N-l)
    den=den*1/N
            
    return num/den

def create_sample(samples, par, ell):
    '''
    Given the pandas dictionary samples, a parameter (such as Ctt) and a multipole number ell, this function returns a list with all samples of Ctt for that multipole number thoughout all samples
    '''

    Cl=samples[par]
    sample=[]
    for i in range(len(samples["Sample #"])):
        if samples['l'][i]==ell:
            sample.append(Cl[i])

    print("\nDebug create_sample:")
    print("# of samples=", max(samples["Sample #"])+1)
    print("len(sample)=", len(sample))

    return sample 

def save_plot(ls, ac, ell, which):
    plt.figure()
    plt.plot(ls, ac)
    plt.title(r"Autocorrelation for {0} at $\ell={1}$".format(which, ell))
    plt.xlabel("l")
    plt.ylabel("Autocorrelation")
    plt.savefig("ac_"+which+"_ell{}.png".format(ell))
    plt.close()


def main():
    sample_file=input('Sample file (with path if necessary): ')
    samples=pd.read_csv(sample_file, sep=' ', header=None, index_col=False, names=['Sample #', 'l', 'Ctt', 'Ctg', 'Cgt', 'Cgg']) #Not sure if this is the correct order
    lmax=int(input("Max. value of lag: "))

    ell_list=[2,3,4,10,20,30,40,50]
    for ell in ell_list:
        Ctt, Cgg, Ctg=create_sample(samples, "Ctt", ell), create_sample(samples, "Cgg", ell), create_sample(samples, "Ctg", ell)    
        ac_tt, ac_gg, ac_tg=[],[],[]
        ls=[]

        for l in range(lmax): 
            ls.append(l)
            ac_tt.append(ac(l, Ctt))
            ac_gg.append(ac(l, Cgg))
            ac_tg.append(ac(l, Ctg))

        save_plot(ls, ac_tt, ell, "Cltt")
        save_plot(ls, ac_gg, ell, "Clgg")
        save_plot(ls, ac_tg, ell, "Cltg")

if __name__=='__main__':
    main()

