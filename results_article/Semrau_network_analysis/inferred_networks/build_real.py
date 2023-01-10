import numpy as np

# We build the real network from the edges for which ChIp-sec interactions
# have been found in the litterature. This concerns only the Stimulus (RA)
# and the three genes Pou5f1, Sox2 and Jarid2, for which such data are available.

p = -np.ones((42, 42))
p[0, 1:] = 0
p[0, [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 20, 22, 23, 25, 27, 29, 34, 35, 36, 37]] = 1

p[1, 1:] = 0
p[1, [1, 2, 3, 4, 5, 6, 7, 10, 12, 20, 22, 27, 33, 37]] = 1

p[2, 1:] = 0
p[2, [1, 2, 3, 4, 7, 9, 10, 13, 21, 25, 35]] = 1

p[6, 1:] = 0
p[6, [3, 7, 10, 11] + list(np.arange(15, 34)) + list(np.arange(35, 42))] = 1


np.save(f'REAL', p)
