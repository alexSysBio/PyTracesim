import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import warnings
import warnings
warnings.filterwarnings('ignore')



class unbound_particle_fast_simulation(object):
    """
    This code is used to simulate the particle displacement in a cell with a nucleoid.

    
    This class includes the following functions:
        show_particle: shows the particle trajectory and the cell/nucleoid mesges if specified
    """



    def __init__(self, interval=0.05, frames=300):
        """
        Initializes a cell and creates a cell polygon.
        
        Input:
            particle_radius: non-negative slope (if zero the particle radius is not considered in the simulation)
            interval: the time interval of the simulation in seconds
        Returns:
            self.dx_list : the list of Δx displacements
            self.dy_list : the list of Δy displacements
            self.particle_radius : the specified particle radius
            self.added_cell: 1 indicating that no cell was added
            self.added_nucleoid : 0 indicating that no nucleoid was added
            self.positions_x_list: the list of the particle positions in the x dimension - used for plotting
            self.positions_y_list: the list of the particle positions in the y dimension - used for plotting
        """

        self.interval = interval
        self.frames = frames
#         self.diffusion_coefficient = diffusion_coefficient
        
    
    def get_particle_displacements(self, diffusion_coefficient):
        
        std_normal = np.sqrt(2*diffusion_coefficient*self.interval)
        dx_random = np.random.normal(0, std_normal, self.frames-1)
        dy_random = np.random.normal(0, std_normal, self.frames-1)
        
        self.dx_array = dx_random
        self.dy_array = dy_random
        
   
    def get_displacement_arrays(self):
        
        return self.dx_array, self.dy_array
    
    
    def get_particle_positions(self):
        
        x_positions = np.cumsum(self.dx_array)
        y_positions = np.cumsum(self.dy_array)
        
        self.positions_x_array = np.append(0, x_positions)
        self.positions_y_array = np.append(0, y_positions)
    
    
    def plot_particle_trajectory(self):
        
        x = list(self.positions_x_array)
        y = list(self.positions_y_array)
        f = list(range(self.frames))
#         print(len(f), len(x), len(y))
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(f), clip=True)
        mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        traj_color = np.array([(mapper.to_rgba(v)) for v in range(len(f))])
        
        for fr in range(self.frames):
            plt.plot(x[fr:fr+2], y[fr:fr+2], color=traj_color[fr])
        
        plt.xlabel('X coordinates (px)')
        plt.ylabel('Y coordinates (px)')
        plt.show()

    
    def run_simulation(self, number_of_particles, diffusion_coefficient_array, seed_number, show):
        
        np.random.seed(seed_number)
        
        mean_dc = np.mean(np.log10(diffusion_coefficient_array))
        std_dc = np.std(np.log10(diffusion_coefficient_array))
        
        print('generating trajectories from a normal distribution of effective')
        print('diffusion coefficients with a mean of', 10**mean_dc,'μm2/sec')
        
        sampled_diffusion_coefficients = 10**np.random.normal(mean_dc, std_dc, number_of_particles)
        
        i = 0
        for par in range(number_of_particles):
            self.get_particle_displacements(sampled_diffusion_coefficients[i])
            self.get_particle_positions()
            
            if show == True:
                print('Diffusion coefficient:',sampled_diffusion_coefficients[i],'μm2/sec')
                self.plot_particle_trajectory()
            
            x = list(self.positions_x_array)
            y = list(self.positions_y_array)
            par_id = 'particle_'+str(par+1)
            
            if i == 0:
                final_df = pd.DataFrame()
                final_df['x'] = x
                final_df['y'] = y
                final_df['par_id'] = par_id
                final_df['par_int'] = i+1
                final_df['time'] = np.arange(0,len(x))*self.interval
                final_df['frame'] = np.arange(0,len(x))
                final_df['dc_sampled'] = sampled_diffusion_coefficients[i]
            elif i >0:
                pre_df= pd.DataFrame()
                pre_df['x'] = x
                pre_df['y'] = y
                pre_df['par_id'] = par_id
                pre_df['par_int'] = i+1
                pre_df['time'] = np.arange(0,len(x))*self.interval
                pre_df['frame'] = np.arange(0,len(x))
                pre_df['dc_sampled'] = sampled_diffusion_coefficients[i]
                final_df = pd.concat([final_df, pre_df])
            
            i+=1
            if i%1000 == 0:
                print(i, 'out of', number_of_particles, 'trajectories')
            
        return final_df
            