import numpy as np

# Setting simulation parameters. This should be the same as in simulation_caller.py

# # Setting simulation parameters
l = np.array([1.0])
td_std = np.array([0.0])  # 1
lambda_std = np.array([0.1])  # 1
delta = np.array([1.0])
beta = np.linspace(0.025, 0.5, num=20)  # 20
r = beta/(1-beta)
alpha = np.array([1.0])  # 1
num_rep = 100  # number of repeats 100
epsilon = np.linspace(0.0, 0.02, 6)  # 6
exp_n = np.array([2,4]) # 2
# epsilon=np.array([0.0]) # no penalty for now.
# # should give 24000 repeats. Run with 1200 job array.
# # should take around 4 hours (240min).

par_vals = {'nstep': 900, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 13, 'delta': delta[0], 'lambda': l[0], 'lambda_std': lambda_std[0], 'td_std':td_std[0]}
X = [len(epsilon), len(exp_n), len(beta), len(alpha), num_rep]
temp1 = np.arange(np.prod(X))  # total number of simulations that need doing
output=np.random.permutation(temp1)
np.save('/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200328_revision_penalty_model/outputs/derangements_1.npy',output)
