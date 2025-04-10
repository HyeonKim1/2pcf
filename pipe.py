import h5py 
import numpy as np
from tqdm import tqdm
from mpi4py import MPI

comm = MPI.COMM_WORLD


class twopcf:
    def __init__(self,data,random_scale=20,boxscale=None):
        self.data=data
        self.N_data=len(data)
        self.N_random=random_scale*len(data)
        if boxscale==None:
            self.min_scale=int(np.min(data))
            self.max_scale=int(np.max(data)+1)
        else:
            self.min_scale=boxscale[0]
            self.max_scale=boxscale[1]
        self.random_data=np.random.uniform(self.min_scale,self.max_scale,size=(int(random_scale*self.N_data),3))
        self.hist_range=np.array([self.min_scale,self.max_scale])

    def Natural(self,bins=10):
        self.DD,array=self.cal_corr_same_data(self.data,bins,self.hist_range)
        print('make DD done')
        self.RR,array=self.cal_corr_same_data(self.random_data,bins,self.hist_range)
        print('make RR done')

        return (self.DD/self.RR)-1 ,array
    
    def DP(self,bins=10):
        self.DD,array=self.cal_corr_same_data(self.data,bins,self.hist_range)
        print('make DD done')
        self.DR,array=self.cal_corr_DR(self.data,self.random_data,bins,self.hist_range)
        print('make DR done')
        return (self.DD/self.DR)-1, array
    
    def Hamilton(self,bins=10):
        
        self.DD,array=self.cal_corr_same_data(self.data,bins,self.hist_range)
        print('make DD done')
        self.RR,array=self.cal_corr_same_data(self.random_data,bins,self.hist_range)
        print('make RR done')
        self.DR,array=self.cal_corr_DR(self.data,self.random_data,bins,self.hist_range)
        print('make DR done')

        return (self.DD*self.RR/(self.DR**2))-1, array
    
    def LS(self,bins=10):
        
        self.DD,array=self.cal_corr_same_data(self.data,bins,self.hist_range)
        print('make DD done')
    
        self.RR,array=self.cal_corr_same_data(self.random_data,bins,self.hist_range)
        print('make RR done')
    
        self.DR,array=self.cal_corr_DR(self.data,self.random_data,bins,self.hist_range)
        print('make DR done')

        return (self.DD-2*self.DR+self.RR)/self.RR, array
    
    def cal_corr_same_data(self,data,bins,hist_range):
        rank = comm.Get_rank()
        size = comm.Get_size()

        local_size = int(len(data) / size)
        extra = np.mod(len(data), size)

        if (rank < extra):
            local_size = local_size + 1
            start_idx = rank * local_size
        else:
            start_idx = rank * local_size + extra 

        end_idx = start_idx + local_size - 1

        hist=np.zeros(bins)
        if rank == 0:
            iterator = tqdm(range(start_idx, end_idx), desc="Progress")
        else:
            iterator = range(start_idx, end_idx)

        for ii in iterator:
            distance=np.sqrt(np.sum((data[ii:]-data[ii])**2,axis=1))
            hist1,array=np.histogram(distance,bins=bins,range=hist_range)
            hist = hist+hist1
        array = (array[:-1] + array[1:]) / 2
        global_hist=comm.allreduce(hist, op=MPI.SUM)

        return global_hist/np.sum(global_hist), array
    
    def cal_corr_DR(self,data,random_data,bins,hist_range):
        rank = comm.Get_rank()
        size = comm.Get_size()

        local_size = int(len(data) / size)
        extra = np.mod(len(data), size)

        if (rank < extra):
            local_size = local_size + 1
            start_idx = rank * local_size
        else:
            start_idx = rank * local_size + extra 

        end_idx = start_idx + local_size - 1

        hist=np.zeros(bins)
        if rank == 0:
            iterator = tqdm(range(start_idx, end_idx), desc="Progress")
        else:
            iterator = range(start_idx, end_idx)

        for ii in iterator:
            distance=np.sqrt(np.sum((random_data-data[ii])**2,axis=1))
            hist1,array=np.histogram(distance,bins=bins,range=hist_range)
            hist = hist+hist1
        
        global_hist=comm.allreduce(hist, op=MPI.SUM)
        array = (array[:-1] + array[1:]) / 2

        return global_hist/np.sum(global_hist), array

        
class read_Gadget4:
    def __init__(self,filename):
        self.filename = filename
        self.my_snapshot=h5py.File(filename)
        self.boxsize=self.my_snapshot['Header'].attrs['BoxSize']
        self.redshift=self.my_snapshot['Header'].attrs['Redshift']

    def read_pos(self,particle_name='CDM'):
        if particle_name=='CDM':
            particle_type='PartType1'
            particle_type_num=1

        else: 
            particle_type='PartType0'
            particle_type_num=0

        self.particle_num=self.my_snapshot['Header'].attrs['NumPart_Total'][particle_type_num]
        self.pos=self.my_snapshot[particle_type]['Coordinates'][:,:]
        return self.pos