import os
from shutil import copyfile
import numpy as np

job_number = 50150767
temp_base_path = '/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200329_revision_correlations_variable_lambda/'
if not os.path.exists(temp_base_path+'outputs_backup'):
    os.mkdir(temp_base_path+'outputs_backup')
# this allows us to know which ones we are re-doing
f = open(temp_base_path+'/failed_{0}.txt'.format(job_number), 'r') # where file_object is the variable to add the file object.
val=[]
for line in f:
    val.append(np.int(line))
# val gives the array number of each one that must be repeated from scratch.
for num in val:
    # first we make a backup copy of these files just in case they are useful
    copyfile(temp_base_path+'outputs/output_vals_{0}.npy'.format(num),
             temp_base_path+'outputs_backup/output_vals_{0}.npy'.format(num))
    copyfile(temp_base_path + 'outputs/progress_tracker_{0}.npy'.format(num),
             temp_base_path + 'outputs_backup/progress_tracker_{0}.npy'.format(num))
    # next we remove these files so we can start from scratch
    os.remove(temp_base_path+'outputs/output_vals_{0}.npy'.format(num))
    os.remove(temp_base_path + 'outputs/progress_tracker_{0}.npy'.format(num))