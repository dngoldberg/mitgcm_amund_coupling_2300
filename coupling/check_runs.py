import io
import subprocess
import pandas as pd

result = subprocess.run('squeue -u dngoldbe -o "%.10i %.20j %.8T"',shell=True, text=True, capture_output=True)
df = pd.read_csv(io.StringIO(result.stdout), delim_whitespace=True)

inds = df['NAME'].str.contains('80')
df = df[inds]
inds = df['STATE'] == "RUNNING"
df = df[inds]
jobs = df.JOBID.values
names = df.NAME.values

print(df['NAME'])

for i in range(len(jobs)):
 subprocess.run('tail -n 2 scripts/slurm-' + str(jobs[i]) + '.out',shell=True, text=True)
 print('run_oce_2009_' + names[i][1:-10] + '_80_iceParm3_coul')
 subprocess.run('tail run_oce_2009_' + names[i][1:-10] + '_80_iceParm3_coul/STDOUT.0001 -n 6',shell=True, text=True)

