import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import growth_simulations as g
import time
import os
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
import seaborn as sns
from sklearn.neighbors import KernelDensity
import scipy.integrate as integrate

temp_ind = np.int(os.environ['SLURM_ARRAY_TASK_ID'])  # extract the ID of this job number and use that to track where we start
temp_ind1 = np.int(os.environ['SLURM_ARRAY_TASK_MAX']) # extract the total number

# temp_ind = 1
print temp_ind, temp_ind1


l = np.array([1.0])
td_std = np.array([0.0])
l_std = np.array([0.2])
delta = np.array([1.0])
beta = np.linspace(0.05, 0.5, num=10) # 9
r = beta/(1-beta)
alpha = np.linspace(0.5, 1.0, num=6) # 5
tic=time.clock()
num_rep = 100
# 8583616  . Run with 1000 job array, 6 repeats per unit.

X = [len(beta), len(alpha), num_rep]
num_celltypes = 3  #mother, daughter, full_pop
a_shape = X+[2]  # zero for simulations, one for numerical estimation based on E-L equation.
print a_shape
celltype = ['Mother', 'Daughter', 'Population']
params = {'nstep': 1200, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 7, 'delta': delta[0], 'lambda': l[0],
            'td_std': td_std[0], 'lambda_std': l_std[0]}

bw = 0.001
L=np.linspace(0.965,1.0,141)
tic=time.time()

temp1 = np.prod(X)  # total number of simulations that need doing
temp2 = temp1/temp_ind1  # number of simulations that will be done by each node
temp_base_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_numerical_estimation_euler_lotka/outputs/'
temp_inds = np.load(temp_base_path+'derangements_1.npy')  # this randomizes which machine gets which simulations and avoids individual ones taking dramatically longer than everyone else (since some asymmetry levels take longer than others for example)
if len(temp_inds)!=temp1:
	raise ValueError('You need to seed me with the right randomized number of simulations')
print 'Must do {0} simulations out of {1}'.format(temp2, temp1)
print X, td_std, l_std, beta, alpha, num_rep
num_celltypes = 3  # mother, daughter, full_pop
observables = ['Sim GR vol', 'EL GR']
a_shape = X+[len(observables)]  # this is the shape of the array that we are appending our results to

if temp_ind<temp_ind1:  # if this is a "regular" job
    temp3 = range((temp_ind-1)*temp2,temp_ind*temp2)
elif temp_ind == temp_ind1:  # if this is the final job in the array
    temp3 = range((temp_ind-1)*temp2,temp1)

prog_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_numerical_estimation_euler_lotka/outputs/progress_tracker_{0}.npy'.format(temp_ind)
if os.path.exists(prog_path):
    progress = np.load(prog_path)[0]
    temp3 = range(progress+1, max(temp3)+1)  # if we have already made some progress in this we must start from where we
    # left off
out_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_numerical_estimation_euler_lotka/outputs/output_vals_{0}.npy'.format(temp_ind)
if os.path.exists(out_path):
    output_val = np.load(out_path)
else:
    output_val = np.zeros(a_shape)
tic=time.time()
print temp3


for i0 in temp3:
    inds=np.unravel_index(i0,X)
    # setting the parameters
    params['r']=beta[inds[0]] / (1 - beta[inds[0]])
    params['alpha']=alpha[inds[1]]
    # doing the simulation
    init_pop = g.starting_popn(params)

    params['nstep'] = 500  # seeding the population with a simulated one
    c, obs, [temp_vols, temp_vols_G1] = g.discr_time_1(params, init_pop)
    init_pop = g.starting_popn_seeded(c, params)
    params['nstep'] = 900  # now we run this simulation for longer with a better seeded population.
    c, obs, [temp_vols, temp_vols_G1] = g.discr_time_1(params, init_pop)

    # calculating the simulation population growth rate and storing it
    temp = scipy.stats.linregress(obs[1][400:], np.log(obs[12][400:]))
    output_val[inds[0],inds[1],inds[2],0]=temp[0]
    # calculating the kernel density estimation
    kde1= KernelDensity(bandwidth=bw,kernel='gaussian')
    kde2= KernelDensity(bandwidth=bw,kernel='gaussian')
    t1=np.array([obj.t_grow for obj in c if obj.celltype==0 and obj.tb>400* params['dt'] * np.log(2)/params['lambda']]).reshape(-1,1)
    t2=np.array([obj.t_grow for obj in c if obj.celltype==1 and obj.tb>400* params['dt'] * np.log(2)/params['lambda']]).reshape(-1,1)
    kde1.fit(t1)
    kde2.fit(t2)
    y1=np.zeros(len(L))
    y2=np.zeros(len(L))
    t1m,t1std=np.mean(t1),np.std(t1)
    xv1=np.concatenate((np.linspace(np.amax([0.0,t1m-3*t1std]),t1m+10*t1std, 5001), np.array(t1m+10*t1std*np.exp(np.linspace(np.log(5)/1000,np.log(5),1000)))),axis=0)
    t2m,t2std=np.mean(t2),np.std(t2)
    xv2=np.concatenate((np.linspace(np.amax([0.0,t2m-3*t2std]),t2m+10*t2std, 5001), np.array(t2m+10*t2std*np.exp(np.linspace(np.log(5)/1000,np.log(5),1000)))),axis=0)
    for temp in range(len(L)):
        y1[temp]=integrate.simps(np.exp(kde1.score_samples(xv1[:,None]))*np.exp(-L[temp]*xv1),xv1)
        y2[temp]=integrate.simps(np.exp(kde2.score_samples(xv2[:,None]))*np.exp(-L[temp]*xv2),xv2)
        if np.mod(temp,10)==0:
            print 'Done L:', temp
    output_val[inds[0],inds[1],inds[2],1]=L[np.argmin(np.absolute(y1/(1-y2)-1))] # calculating the value of L that fits
    # for these parameter values
    print 'I have done {0} parameter sets.'.format(i0+1), 'time taken:',time.time()-tic
    print 'beta',beta[inds[0]], 'alpha',alpha[inds[1]], 'repeat {0}'.format(inds[2])
    print output_val[inds[0],inds[1],inds[2],:]

    # saving the results so that we can restart from where we left off
    np.save(out_path, output_val)
    # showing that we have made this progress
    progress = np.array([i0])
    np.save(prog_path, progress)
    if np.mod(i0 - temp3[0], 10) == 0:
        print 'I have completed: {0} repeats. Time taken ='.format(i0 - temp3[0] + 1), time.time() - tic
    del c, obs, init_pop, t1, t2
print 'I have completed: {0} (all) repeats. Time taken ='.format(len(temp3)), time.time() - tic