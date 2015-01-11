"""Module for reading and

The data for this comes from the Saguaro Astronomy Club Database,
which contains many deep sky objects along with statistics useful to
any amateur astronomer."""

import pdb, os.path
import body

inputfilename = 'SAC_DeepSky_770_Fence.txt'
inputpath = os.path.join(os.path.dirname(__file__), inputfilename)
inputfile = file(inputpath)

pdb.set_trace()
# Skip over the header line.
for line in inputfile:
    if line.strip():
        break
# Process the regular lines.
for line in inputfile:
    line = line.strip()
    if not line: continue
    
    
