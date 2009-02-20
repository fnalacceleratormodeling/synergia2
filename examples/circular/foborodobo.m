# calculations going into fobrdobo_s.lat for rf cavity

format("long")

pct = 0.4

momentum = 100 # gev/c
bg = momentum/mass_proton # beta*gamma lorentz factor
g = sqrt(bg**2 + 1) # gamma lorentz factor
b = bg/g # beta lorentz factor

n = 8 # cells in the ring
bendangle = 2*pi/n
focus = 7 # focal length of quad
sepn = 10 # distance between quad centers
quadlength = 0.2 # length of quad
strength = 1/(focus*quadlength)
pct = 0.4 # fraction of space between quads occupied by dipoles
bendlength = pct * (sepn - quadlength)
driftlength = (sepn - quadlength - bendlength)/2
harmno = 32 # harmonic number should be multiple of n
				# in order for it to do a full cycle
				# when a particle gets to the next rf
				# cavity.
cell_length = quadlength + driftlength + bendlength + \
    driftlength + quadlength + driftlength + bendlength + driftlength

rf_freq = harmno * speed_of_light * b/(n * cell_length)
