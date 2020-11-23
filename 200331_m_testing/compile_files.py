import numpy as np

# num_array = 2200
num_array = 1000
path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200331_m_testing/outputs/'
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

# scp fbarber@login.rc.fas.harvard.edu:/n/regal/murray_lab/fjbarber/19_growth_rate_project/200323_revision_no_penalty_model/outputs/output_compiled.npy ~/Documents/18_growth_rate_project/cluster_running/

# scp ./200323_revision_no_penalty_model/* fbarber@login.rc.fas.harvard.edu:/n/murraylab/Users/fbarber/20_growth_rate_project/200323_revision_no_penalty_model/