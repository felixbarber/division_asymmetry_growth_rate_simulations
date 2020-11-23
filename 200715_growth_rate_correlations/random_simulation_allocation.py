import numpy as np

# Setting simulation parameters. This should be the same as in simulation_caller.py

# # Setting simulation parameters
l = np.array([1.0])
td_std = np.linspace(0.0, 0.2, 3)  # 3
lambda_std = np.linspace(0.0, 0.2, num=3)  # 3
delta = np.array([1.0])
beta = np.linspace(0.1, 0.5, num=17)  # 17
r = beta/(1-beta)
alpha = np.linspace(0.5, 1.0, num=2)  # 2
num_rep = 50  # number of repeats
gr_corr_a=np.linspace(0.0,0.9,num=10)  # 10
# # should give 153000 repeats. Run with 1500 job array.
# # should take around 15 hours. 65315013

par_vals = {'nstep': 900, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 16, 'delta': delta[0], 'lambda': l[0]}
X = [len(td_std), len(lambda_std), len(beta), len(alpha), num_rep, len(gr_corr_a)]
temp1 = np.arange(np.prod(X))  # total number of simulations that need doing
output=np.random.permutation(temp1)
np.save('/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200715_growth_rate_correlations/outputs/derangements_1.npy',output)
