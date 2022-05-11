import numpy as np
import scipy.io as io
import xarray as xr

class micom_tools:
    def __init__(self):
        """
        Initialize the reader with fixed grid size.
        """
        # Define size
        self.Lx = 1200   # km
        self.Ly = 191.25 # km
        self.D = 4300    # m
        self.scalingm = 98060
        self.scalingcm = 980.6

        # Indices
        self.jdd=322
        self.idd=53
        self.k1=43
       
        self.ii=321
        self.jj=52
        self.kdd=42
        self.kk=44
        
        self.dx = self.Lx / self.jdd
        self.dy = self.Ly / self.idd
        self.dz = self.D / self.k1
        
        # Arrays with lengths
        self.C = np.arange(self.jdd) * self.dx + self.dx/2
        self.yc = np.arange(self.idd) * self.dy + self.dy/2
        self.zc = np.arange(self.k1) * self.dz + self.dz/2
        
        self.xp = np.arange(self.jdd + 1) * self.dx
        self.yp = np.arange(self.idd + 1) * self.dy
        self.zp = np.arange(self.k1 + 1) * self.dz

    def read_output(self, path, returntype='xarray'):
        """
        Read binary output files from the MICOM model.
        
        Arguments
        ---------
        path : str
            Path to file
        
        returntype : str
            `xarray` or `numpy`
        
        Returns
        -------
        np.array
            3D array with shape (k1, jdd, idd) = (43, 322, 53)
        or
        xr.DataArray
            3D DataArray with basin dimensions
            
        """
        f = io.FortranFile(path, "r")
        output = np.flip(f.read_reals(np.dtype('<f8')).reshape(self.k1,self.jdd,self.idd), 0)
        output = np.flip(output, 2)
        f.close()
        
        if returntype == 'numpy':
            return output
        elif returntype == 'xarray':
            output = xr.DataArray(data = output, 
                                  coords = [self.zc, self.C, self.yc],
                                  dims = ['Z', 'X', 'Y'])
            output['X'].attrs = {"long_name" : "X-coordinate of the cell center",
                                 "unit" : "kilometers"}
            output['Y'].attrs = {"long_name" : "Y-coordinate of the cell center",
                                 "unit" : "kilometers"}
            output['Z'].attrs = {"long_name" : "Z-coordinate of the cell center",
                                 "unit" : "meters"}
            return output
        else:
            raise ValueError("`returntype` must be either `numpy` or `xarray`.")
    
    
    def merge_along_time(self, ds_list, dt=1.5):
        """
        Merge xarray DataArrays along a time dimension.
        
        Arguments
        ---------
        ds_list : list, tuple
            Iterable of xr.DataArray items, output by `read_output()`
        dt : int, float
            timestep between files in hours
        """
        ds = xr.concat(ds_list, dim='time')
        n_time = ds.time.size
        ds = ds.assign_coords(time=np.arange(0, n_time*dt, dt))
        ds['time'].attrs = {"long_name" : "time referenced to 0",
                            "unit" : "hours"}
        return ds
        
    
    def harmonic_analysis(self, analyses, average, adjust_phase=False):
        """
        Compute harmonic analysis using analysis files and an average file.
        
        Arguments
        ---------
        analyses : xr.DataArray
            DataArray that contains several timesteps
        average : xr.DataArray
            DataArray that contains the average of the output
        adjust_phase : bool
            Shift the phase by pi. This may be necessary to get to a similar phase as in the paper.
        """
        
        n = analyses.time.size
        dt = analyses.time[1] - analyses.time[0] # should be in hours
        anomalies = (analyses - average) / self.scalingcm
        acum = np.cumsum(anomalies, axis=1)
        mean = acum.mean('time')
        sumsq = (acum**2).sum('time')
        residual = 0. * mean
        C = 0. * mean
        S = 0. * mean
        for l in range(0, n):
            # Here we generalize the computation to work with any dt
            C += (acum[l, :, :, :] - mean) * np.cos((l+1)*np.pi / (24/dt/4))
            S += (acum[l, :, :, :] - mean) * np.sin((l+1)*np.pi / (24/dt/4))
        C = C / (n/2)
        S = S / (n/2)
        for l in range(0, n):
            residual += (mean + C * np.cos((l+1) * np.pi/4) + S * np.sin((l+1) * np.pi/4) - acum[l, :, :, :])**2
        residual = (residual/sumsq).drop('time')
        amplitude = (np.sqrt(C**2 + S**2)).drop('time') / self.scalingcm
        phase = (np.arctan2(S, C)).drop('time')
        
        residual.name = 'Residual [cm]'
        amplitude.name = 'Amplitude [cm]'
        phase.name = 'Phase [rad]'
        
        if adjust_phase:
            phase = phase + np.pi - ((phase > 0).astype('float')*np.pi*2)
        
        return amplitude, phase, residual
