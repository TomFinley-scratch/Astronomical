import sys, pdb

# Ignore stars with visual magnitude greater than 9.
cutoff = 9

hip_data_path = '/Volumes/Pallas/Hipparcos Star Data/hip_main.dat'

f = file('/Volumes/Pallas/Tycho-2/catalog.dat')
processed = 0
min_mag = 20
hip_identified = 0
for line in f:
    tokens = line.split('|')
    try:
        hip_number = int(tokens[23])
        hip_identified += 1
    except ValueError:
        pass
    try:
        visual_magnitude = float(tokens[19])
    except ValueError:
        continue
    if visual_magnitude < min_mag:
        print
        print 'Found min', visual_magnitude
        min_mag = visual_magnitude
        if visual_magnitude < 3:
            print line
    processed += 1
    if (processed%1000)==0:
        sys.stdout.write('.')
        sys.stdout.flush()
f.close()
print
print 'Hip identified:', hip_identified
