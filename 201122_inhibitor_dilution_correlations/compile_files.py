import numpy as np

# num_array = 2200
num_array = 2400
path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_inhibitor_dilution_correlations/outputs/'
for i0 in range(1,num_array+1):
    name = 'output_vals_{0}.npy'.format(i0)
    temp=np.load(path+name)
    if i0==1:
        output=temp
    else:
        output+=temp
    if np.mod(i0,100)==0:
        print "I have done {0} output arrays".format(i0)
out_name='output_compiled.npy'
np.save(path+out_name,output)

for i0 in range(1,num_array+1):
    name = 'progress_tracker_{0}.npy'.format(i0)
    temp=np.load(path+name)
    if i0==1:
        output=np.zeros(num_array)
    output[i0-1]=temp[0]
    if np.mod(i0,100)==0:
        print "I have done {0} progress arrays".format(i0)
out_name1='progress_compiled.npy'
np.save(path+out_name1,output)

# scp fbarber@login.rc.fas.harvard.edu:/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201122_inhibitor_dilution_correlations/outputs/output_compiled.npy ~/Documents/18_growth_rate_project/cluster_running/

# scp ./201007_correlations_validation_theory_model_budding/* fbarber@login.rc.fas.harvard.edu:/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201007_correlations_validation_theory_model_budding/