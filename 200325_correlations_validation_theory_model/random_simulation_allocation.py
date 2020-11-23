import numpy as np

# Setting simulation parameters. This should be the same as in simulation_caller.py

# # Setting simulation parameters
l = np.array([1.0])
td_std = np.linspace(0.0, 0.1, 5)  # 3
lambda_std = np.linspace(0.0, 0.1, num=2)  # 3
delta = np.array([1.0])
beta = np.linspace(0.025, 0.5, num=20)  # 20
r = beta/(1-beta)
alpha = np.linspace(0.2, 1.0, num=5)  # 11
num_rep = 50  # number of repeats
# # should give 50000 repeats. Run with 1250 job array.
# # should take around 7 hours.


par_vals = {'nstep': 1000, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 7, 'delta': delta[0], 'lambda': l[0]}
X = [len(td_std), len(lambda_std), len(beta), len(alpha), num_rep]
temp1 = np.arange(np.prod(X))  # total number of simulations that need doing
output=np.random.permutation(temp1)
np.save('/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200325_correlations_validation_theory_model/outputs/derangements_1.npy',output)
