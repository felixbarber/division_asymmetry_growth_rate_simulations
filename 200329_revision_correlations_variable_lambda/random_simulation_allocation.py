import numpy as np

# Setting simulation parameters. This should be the same as in simulation_caller.py

# # Setting simulation parameters
l = np.array([1.0])
td_std = np.linspace(0.0, 0.3, 4)  # 4
lambda_std = np.linspace(0.0, 0.3, num=7)  # 7
delta = np.array([1.0])
beta = np.linspace(0.025, 0.5, num=20)  # 20
r = beta/(1-beta)
alpha = np.array([0.5])  # 1
num_rep = 100  # number of repeats
# # should give 56000 repeats. Run with 1400 job array.
# # should take around 7 hours. 600min.

par_vals = {'nstep': 900, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 15, 'delta': delta[0], 'lambda': l[0]}
X = [len(td_std), len(lambda_std), len(beta), len(alpha), num_rep]
temp1 = np.arange(np.prod(X))  # total number of simulations that need doing
output=np.random.permutation(temp1)
np.save('/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200329_revision_correlations_variable_lambda/outputs/derangements_1.npy',output)