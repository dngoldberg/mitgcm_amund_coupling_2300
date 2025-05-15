from mds import rdmds, parsemeta
import numpy as np
import re
import os
import sys

def writemeta(data, filename):
    """Writes a dictionary to a file in the format expected by parsemeta()."""
    with open(filename, 'w') as f:
        for key, values in data.items():
            if all(isinstance(v, str) for v in values):
                # Write string list in {} format
                values_str = ' '.join(f"'{v}'" for v in values)
                f.write(f"{key} = {{{values_str}}};\n")
            else:
                # Write numeric or mixed list in [] format
                values_str = ', '.join(str(v) for v in values)
                f.write(f"{key} = [{values_str}];\n")

def average_timesteps (name, coupling_period, averaging_period, final_time, dt): 

    # name: prefix name
    # coupling_period: in timesteps, the expected period of the diag file
    # averaging_period: in timesteps, the period over which we wish to average
    # final_time: timestep of last file
    # dt: timestep in sec

    number_of_files = averaging_period/coupling_period

    if (number_of_files != np.round(number_of_files)): 
        raise ValueError('averaging period should be a multiple of coupling period')

    nfiles = int(number_of_files)

    arrs, x, meta = rdmds(name, final_time, returnmeta=True)
    metafilename = name + '.' + str(final_time).zfill(10) + '.meta'

    meta2 = parsemeta(metafilename)
    time2 = meta2['timeInterval'][0]
    meta2['timeInterval'] = [time2-averaging_period*dt, time2]
    writemeta(meta2, metafilename)

    arrnew = arrs / number_of_files

    for i in np.arange(final_time-coupling_period, final_time-averaging_period, -coupling_period):
        
        arrtime = rdmds(name, i, returnmeta=True)
        arrs = arrs + arrtime / number_of_files
        os.remove(name + '.' + str(i).zfill(10) + '.meta')
        os.remove(name + '.' + str(i).zfill(10) + '.data')

    arrs.byteswap().tofile(name + '.' + str(final_time).zfill(10) + '.data')

    print('averaged diagnostics ' + name + '.' + str(final_time).zfill(10) + '.data')

if __name__ == "__main__":
    name=sys.argv[1]
    coupling_period=int(sys.argv[2])
    averaging_period=int(sys.argv[3])
    final_time=int(sys.argv[4])
    dt = float(sys.argv[5])
    print (sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3] + ' ' + sys.argv[4] + ' ' + sys.argv[5])
    average_timesteps(name, coupling_period, averaging_period, final_time, dt)

