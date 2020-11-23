import numpy as np

# Setting simulation parameters. This should be the same as in simulation_caller.py

l = np.array([1.0])
td_std = np.array([0.0])
l_std = np.array([0.2])
delta = np.array([1.0])
beta = np.linspace(0.05, 0.5, num=10) # 9
r = beta/(1-beta)
alpha = np.linspace(0.5, 1.0, num=6) # 5
num_rep = 100
X = [len(beta), len(alpha), num_rep]
num_celltypes = 3  #mother, daughter, full_pop
a_shape = X+[2]  # zero for simulations, one for numerical estimation based on E-L equation.
print a_shape
celltype = ['Mother', 'Daughter', 'Population']
params = {'nstep': 1200, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 7, 'delta': delta[0], 'lambda': l[0],
            'td_std': td_std[0], 'lambda_std': l_std[0]}

bw = 0.001
L=np.linspace(0.965,1.0,141)

X = [len(beta), len(alpha), num_rep]
temp1 = np.arange(np.prod(X))  # total number of simulations that need doing
output=np.random.permutation(temp1)
np.save('/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_numerical_estimation_euler_lotka/outputs/derangements_1.npy',output)
