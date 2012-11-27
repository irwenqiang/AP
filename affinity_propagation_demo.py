from numpy import *
import numpy as np

damping = 0.5
n_samples=4
max_iter=9
A = np.zeros((n_samples, n_samples))
R = np.zeros((n_samples, n_samples))  

ind = np.arange(n_samples)
S = array([[-0.04,-0.01,-0.04,-0.03], [-0.01, -0.04, -0.02, -0.01], [-0.01, -0.01, -0.04, -0.02], [-0.02, -0.01,-0.02, -0.04]])
for it in range(max_iter):
	# Compute responsibilities
	Rold = R.copy()

	# a(i, k) + s(i, k)
	AS = A + S

	# indices of max(a(i, k) + s(i, k))
	I = np.argmax(AS, axis=1)
	# max(a(i, k) + s(i, k))
	Y = AS[np.arange(n_samples), I]  # np.max(AS, axis=1)

	AS[ind, I[ind]] = - np.finfo(np.double).max
	Y2 = np.max(AS, axis=1)

	# r(i, k) <- s(i,k) - max(a(i, k) + s(i, k))
	R = S - Y[:, np.newaxis]

	# r(i, k) <- s(i,k) - max(a(i, k') + s(i, k'))
	R[ind, I[ind]] = S[ind, I[ind]] - Y2[ind]

	R = (1 - damping) * R + damping * Rold  # Damping
	print "R" + "-" * 50
	print R[0:]

	# Compute availabilities
	Aold = A
	Rp = np.maximum(R, 0)
	Rp.flat[::n_samples + 1] = R.flat[::n_samples + 1]
	#print "Rp:"
	#print Rp[0:]

	A = np.sum(Rp, axis=0)[np.newaxis, :] - Rp
	#print "sum:"
	#print np.sum(Rp, axis=0)[np.newaxis,:]
	#print "A:"
	#print A[0:]
	dA = np.diag(A)
	A = np.minimum(A, 0)

	A.flat[::n_samples + 1] = dA

	A = (1 - damping) * A + damping * Aold  # Damping
	print "A" + "-" * 50
	print A[0:]

