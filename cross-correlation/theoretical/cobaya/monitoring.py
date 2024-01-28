from cobaya.samplers.mcmc import plot_progress

chains_path_name=input("Chains' directory path/chain name (ex: 'my_samples/chain_name'):\n")

plot_progress(chains_path_name)#, fig_args={"figsize": (6,4)})
import matplotlib.pyplot as plt
plt.tight_layout()
plt.show()
