import numpy as np
import growth_simulations as g
import os
import scipy
import time
import weakref


# IMPORTANT NOTE: MAKE SURE YOU DELETE ALL THE CONTENTS OF ./outputs/ BEFORE RUNNING THIS SCRIPT, OTHERWISE IT WILL NOT
# WORK PROPERLY.
temp_ind = np.int(os.environ['SLURM_ARRAY_TASK_ID'])  # extract the ID of this job number and use that to track where we start
temp_ind1 = np.int(os.environ['SLURM_ARRAY_TASK_MAX']) # extract the total number
print temp_ind, temp_ind1

# to select the ones that failed sacct -j 66088472 -o "JobID%16,Partition,State" | grep serial_re | grep -e 'OUT_OF_ME' -e FAILED -e TIMEOUT | awk '{print $1}' | awk -F _ '{print $2}' >> failed_66088472.txt

# in this array. Assumes the arrays start from 1.

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
# # should take around 7 hours. 600min.50150767 

par_vals = {'nstep': 900, 'dt': 0.01, 'v_init': 1.0, 'modeltype': 15, 'delta': delta[0], 'lambda': l[0]}
X = [len(td_std), len(lambda_std), len(beta), len(alpha), num_rep]
temp1 = np.prod(X)  # total number of simulations that need doing
temp2 = temp1/temp_ind1  # number of simulations that will be done by each node
temp_base_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200329_revision_correlations_variable_lambda/outputs/'
temp_inds = np.load(temp_base_path+'derangements_1.npy')  # this randomizes which machine gets which simulations and avoids individual ones taking dramatically longer than everyone else (since some asymmetry levels take longer than others for example)
if len(temp_inds)!=temp1:
	raise ValueError('You need to seed me with the right randomized number of simulations')
print 'Must do {0} simulations out of {1}'.format(temp2, temp1)
print X, td_std, lambda_std, beta, alpha, num_rep
num_celltypes = 3  # mother, daughter, full_pop
observables = ['MD corr', 'MD corr p val', 'GR num', 'GR vol', 'td_std', 't_av', 'vb_std', 'v_av', 'linear slope']
a_shape = X+[num_celltypes, len(observables)]  # this is the shape of the array that we are appending our results to

if temp_ind<temp_ind1:  # if this is a "regular" job
    temp3 = range((temp_ind-1)*temp2,temp_ind*temp2)
elif temp_ind == temp_ind1:  # if this is the final job in the array
    temp3 = range((temp_ind-1)*temp2,temp1)

prog_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200329_revision_correlations_variable_lambda/outputs/progress_tracker_{0}.npy'.format(temp_ind)
if os.path.exists(prog_path):
    progress = np.load(prog_path)[0]
    temp3 = range(progress+1, max(temp3)+1)  # if we have already made some progress in this we must start from where we
    # left off
out_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200329_revision_correlations_variable_lambda/outputs/output_vals_{0}.npy'.format(temp_ind)
if os.path.exists(out_path):
    output_val = np.load(out_path)
else:
    output_val = np.zeros(a_shape)
tic=time.time()
print temp3
for ind in temp3:
    ind1 = np.unravel_index(temp_inds[ind],X)  # this is the index in the saved array that we must go through
    par_vals['td_std'] = td_std[ind1[0]]
    par_vals['lambda_std']=lambda_std[ind1[1]]
    par_vals['r']=r[ind1[2]]
    par_vals['alpha']=alpha[ind1[3]]
    # running the actual simulation
    init_pop = g.starting_popn(par_vals)
    par_vals['nstep']=500  # seeding the population with a simulated one
    c, obs, [temp_vols, temp_vols_G1] = g.discr_time_1(par_vals, init_pop)
    init_pop = g.starting_popn_seeded(c, par_vals)
    par_vals['nstep'] = 900  # now we run this simulation for longer with a better seeded population.
    c, obs, [temp_vols, temp_vols_G1] = g.discr_time_1(par_vals, init_pop)

    # calculating the relevant observables from the simulation results. Note we calculate things after 400 steps based
    # on what was most accurate when varying this burn length (distance from theoretical prediction and spread in
    # repeats).
    for i6 in range(2):
        xv = [obj.t_grow for obj in c if obj.celltype == i6 and obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda'] and obj.parent.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
        yv = [obj.parent.t_grow for obj in c if obj.celltype == i6 and obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda'] and obj.parent.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
        temp = scipy.stats.pearsonr(xv, yv)
        # print temp
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 0] = temp[0]  # PCC
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 1] = temp[1]  # PCC pval
        if i6 == 0:
            temp = scipy.stats.linregress(obs[1][400:], np.log(obs[10][400:]))
            output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 2] = temp[0]  # GR num
            temp = scipy.stats.linregress(obs[1][400:], np.log(obs[12][400:]))
            output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 3] = temp[0]  # GR vol
        if i6 == 1:
            temp = scipy.stats.linregress(obs[1][400:], np.log(obs[11][400:]))
            output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 2] = temp[0]  # GR num
            temp = scipy.stats.linregress(obs[1][400:], np.log(obs[13][400:]))
            output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 3] = temp[0]  # GR vol
        xv = [obj.t_grow for obj in c if obj.celltype == i6 and obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 4] = np.std(xv)  # td_std
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 5] = np.mean(xv)  # td_av
        xv = [obj.vb for obj in c if obj.celltype == i6 and obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 6] = np.std(xv)  # vb_std
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 7] = np.mean(xv)  # vb_av
        yv = [obj.vd for obj in c if obj.celltype == i6 and obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
        temp = scipy.stats.linregress(xv, yv)
        output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 8] = temp[0]  # vb vd linear slope
    i6 = 2
    xv = [obj.t_grow for obj in c if
          obj.celltype == 1 and obj.tb > 400 * par_vals['dt'] * np.log(2) / par_vals['lambda']]  # daughter progeny
    yv = [obj.parent.nextgen.t_grow for obj in c if
          obj.celltype == 1 and obj.tb > 400 * par_vals['dt'] * np.log(2) / par_vals['lambda']]  # mother progeny
    temp = scipy.stats.pearsonr(xv, yv)
    # print temp
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 0] = temp[0]  # PCC
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 1] = temp[1]  # PCC pval

    temp = scipy.stats.linregress(obs[1][400:], np.log(obs[4][400:]))
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 2] = temp[0]  # GR num
    temp = scipy.stats.linregress(obs[1][400:], np.log(obs[7][400:]))
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 3] = temp[0]  # GR num
    xv = [obj.t_grow for obj in c if obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 4] = np.std(xv)  # td_std
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 5] = np.mean(xv)  # td_av
    xv = [obj.vb for obj in c if obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 6] = np.std(xv)  # vb_std
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 7] = np.mean(xv)  # vb_av
    yv = [obj.vd for obj in c if obj.tb>400* par_vals['dt'] * np.log(2)/par_vals['lambda']]
    temp = scipy.stats.linregress(xv, yv)
    output_val[ind1[0], ind1[1], ind1[2], ind1[3], ind1[4], i6, 8] = temp[0]  # vb vd linear slope

    # saving the results so that we can restart from where we left off
    np.save(out_path,output_val)
    # showing that we have made this progress
    progress = np.array([ind])
    np.save(prog_path,progress)
    if np.mod(ind-temp3[0],10)==0:
        print 'I have completed: {0} repeats. Time taken ='.format(ind-temp3[0]+1), time.time()-tic
    del c, obs, init_pop, xv, yv
print 'I have completed: {0} (all) repeats. Time taken ='.format(len(temp3)), time.time()-tic