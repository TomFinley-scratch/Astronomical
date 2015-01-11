constellation_data_path = "constellationship.fab"

def __read():
    """Reads constallation figures by Hipparcos catalog number.

    This will return a list of two element tuples.  The first element
    is a small text identifier for the constallation (e.g., 'Ori' for
    Orion).  The second is a list of duples, with each duple
    containing two Hipparcos catalog numbers indicating that the
    constellation figure has a line between these two.  For example,
    the figure corresponding to Triangulum would be represented by the
    value:

    'Tri', [(10559, 10064), (10064, 8796), (8796, 10559)]

    This indicates that Triangulum is represented as a figure drawn
    between the indicated pairs of stars."""
    f = file(constellation_data_path)
    constellations = []
    for line in f:
        tokens = line.split()
        if not tokens: continue
        hip_numbers = [int(t) for t in tokens[2:]]
        element = tokens[0], zip(hip_numbers[::2], hip_numbers[1::2])
        constellations.append(element)
    f.close()
    return constellations

constellations = __read()
