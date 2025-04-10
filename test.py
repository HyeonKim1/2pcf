from pipe import *
from mpi4py import MPI


filename='./snapshot_002.hdf5'
my_bins=30

my_snapshot=read_Gadget4(filename)
my_pos=my_snapshot.read_pos()

my_2pct=twopcf(my_pos,boxscale=[0,my_snapshot.boxsize])

psi,array=my_2pct.Natural(bins=my_bins)
results={'psi':psi,'array':array}
np.savez('Natural.npz',**results)

DP_psi,DP_array=my_2pct.DP(bins=my_bins)
results={'psi':DP_psi,'array':DP_array}
np.savez('DP.npz',**results)

Hamilton_psi,Hamilton_array=my_2pct.Hamilton(bins=my_bins)
results={'psi':Hamilton_psi,'array':Hamilton_array}
np.savez('Hamilton.npz',**results)

LS_psi,LS_array=my_2pct.LS(bins=my_bins)
results={'psi':LS_psi,'array':LS_array}
np.savez('LS.npz',**results)


