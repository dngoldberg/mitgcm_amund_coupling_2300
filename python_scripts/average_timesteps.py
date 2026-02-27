from mds import rdmds, parsemeta
import numpy as np
import re
import os
import sys
#from IPython import embed

#def writemeta(data, filename):
#    """Writes a dictionary to a file in the format expected by parsemeta()."""
#    with open(filename, 'w') as f:
#        for key, values in data.items():
#            if all(isinstance(v, str) for v in values):
#                # Write string list in {} format
#                values_str = ' '.join(f"'{v}'" for v in values)
#                f.write(f"{key} = {{{values_str}}};\n")
#            else:
#                # Write numeric or mixed list in [] format
#                values_str = ', '.join(str(v) for v in values)
#                f.write(f"{key} = [{values_str}];\n")

def writemeta(data, filepath):
    """Writes dictionary `data` to `filepath` in the custom meta format."""
    with open(filepath, 'w') as f:
        for key, values in data.items():
            if key == 'dataprec':
                # Use square brackets for 'dataprec'
                formatted_values = ', '.join(str(v) for v in values)
                line = f"{key} = ['{formatted_values}'];\n"
            elif all(isinstance(v, str) for v in values):
                # Use curly brackets for lists of strings
                formatted_values = ' '.join(f"'{v}'" for v in values)
                line = f"{key} = {{ {formatted_values} }};\n"
            else:
                # Default: assume it's a list of numbers, use square brackets
                formatted_values = ', '.join(str(v) for v in values)
                line = f"{key} = [{formatted_values}];\n"

            f.write(line)

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
    time2 = meta2['timeInterval'][1]
    meta2['timeInterval'] = [time2-averaging_period*dt, time2]
    print([time2-averaging_period*dt, time2])
    writemeta(meta2, metafilename)

    arrs /= number_of_files

    for i in np.arange(final_time-coupling_period, final_time-averaging_period, -coupling_period):
        
        arrtime = rdmds(name, i, returnmeta=False)
        arrs = arrs + arrtime / number_of_files
        os.remove(name + '.' + str(i).zfill(10) + '.meta')
        os.remove(name + '.' + str(i).zfill(10) + '.data')

    if meta2['dataprec'][0] == 'float32':
       print('compressing to 32 bit')
       newarrs = arrs.astype('float32')
    else:
       print('keeping at 64 bit')
       newarrs = arrs

    newarrs.byteswap().tofile(name + '.' + str(final_time).zfill(10) + '.data')

    print('averaged diagnostics ' + name + '.' + str(final_time).zfill(10) + '.data')

if __name__ == "__main__":
    name=sys.argv[1]
    coupling_period=int(sys.argv[2])
    averaging_period=int(sys.argv[3])
    final_time=int(sys.argv[4])
    dt = float(sys.argv[5])
    print (sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3] + ' ' + sys.argv[4] + ' ' + sys.argv[5])
    average_timesteps(name, coupling_period, averaging_period, final_time, dt)

