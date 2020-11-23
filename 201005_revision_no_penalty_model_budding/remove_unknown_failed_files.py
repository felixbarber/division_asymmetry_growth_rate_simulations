import os
import numpy as np

job_number = 4255716

path1='/n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/201005_revision_no_penalty_model_budding/failed_unknown_reason_{0}.txt'.format(job_number)
# this allows us to know which ones we are re-doing
f = open(path1, 'r') # where file_object is the variable to add the file object.
val=[]
for line in f:
    val.append(np.int(line))
print val

for num in val:
    temp_path ="./outputs/output_vals_{0}.npy".format(num)
    if os.path.exists("demofile.txt"):
        os.remove(temp_path)
    else:
        print("The vals file does not exist at location {0}".format(num))
    temp_path = "./outputs/progress_tracker_{0}.npy".format(num)
    if os.path.exists("./outputs/progress_tracker_{0}.npy".format(num)):
        os.remove(temp_path)
    else:
        print("The progress file does not exist at location {0}".format(num))



