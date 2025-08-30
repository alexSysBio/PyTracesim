import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from shapely.geometry import Point,Polygon,LineString
import random
import warnings
warnings.filterwarnings('ignore')



class bound_particle_simulation(object):
    """
    This code is used to simulate the particle displacement in a cell with a nucleoid.

    
    This class includes the following functions:
        add_cell: adds a cell in the simulation
        add_nucleoid: adds a nucleoid in the cell object
        get_cell_polygon: returns a polygon object of the constructed cell boundaries
        get_nucleoid_polygon: returns a polygon object of the constructed nucleoid boundaries
        get_active_polygon: returns the region in the cell that does not belong into the nucleoid
        is_position_in_cell: checks if a particle position is in the cell or not
        is_position_in_nucleoid: checks if a particle position is in the nucleoid or not
        is_position_in_active_polygon: checks if a particle position is in the active polygon
        nucleoid_crossing: checks if a particle displacement crosses the nucleoid boundaries
        random_point_within_cell: samples one or more random positions from the active polygon
        set_initial_position: sets the initial position of the particle
        get_displacement_list: get the lsit of x/y dosplacements to the most recent time point
        get_position: get the current particle position in the simulation
        set_next_position: set the next particle position in the simulation and update the current position of the particle
        show_particle: shows the particle trajectory and the cell/nucleoid mesges if specified
    """



    def __init__(self, particle_radius=0.1, interval=0.05):
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

        # initialize lists to store the dx and dy displacements
        self.dx_list = []
        self.dy_list = []
        
        self.particle_radius = particle_radius # store the particle radius in the class
        
        self.added_cell = 0  # initialize the class with no cell aded
        self.added_nucleoid = 0 # initialize the class with no nucleoid added
        
        self.positions_x_list = []
        self.positions_y_list = []
        
        self.interval = interval
    
    
    def design_geometrical_shape(self, width, length, show=False):
        """
        This function is used to design a geometrical shape like a bacilus.
        The width and the length of the geometrical shape are specified by the user.
        A semi-circle with a radius half the cell width is used to create the rounded cell ends.
        
        Input:
            width: positive float - the cell width in micrometers
            length: positive float - the cell length in micrometers
        Returns:
            shape_mesh: tuple with two lists, each corresponding to the x and y coordinates of the cell mesh ([x1,x2,..,xn],[y1,y2,...,yn])
            shapely_points: list of lists corresponding to the mesh coordinates, compatible with the shapely library in python ([x1,y1],[x2,y2],...,[xn,yn])
            polygon: a polygon object from the shapely library in python
        
        """
        
        # create empty lists to store the mesh coordinates
        y_coord = []
        x_coord = []
        # create the length with suboixel resolution (0.01 pixels). This allows for a smooth mesh creation
        x_array = np.arange(0,length+0.01,0.01)
#        print(x)
        
        # Iterate across the shape length and generate half the shape along its length
        for x in x_array:
            
            if x < width/2:   # create the semicircle at the left shape edge
                #print(x)
                x_triangle = width/2 - x
                y = np.sqrt(round(((width/2)**2 - x_triangle**2),10))
#                print(y)
                y_coord.append(y)
                x_coord.append(x-length/2)

                
            elif x > length-width/2:  # create a semicircle at the right shape edge
                x_triangle = x - (length-width/2)
                y = np.sqrt(round(((width/2)**2 - x_triangle**2),10))
#                print(y)
                y_coord.append(y)
                x_coord.append(x-length/2)

            else:   # for the rest of the shape length in the middle
                y = width/2
                y_coord.append(y)
                x_coord.append(x-length/2)

        
        # Create the other half of the shape
        # For the X dimension
        x_coord_reversed = x_coord.copy()
        x_coord_reversed = x_coord_reversed[0:-1]  # remove the last point
        x_coord_reversed.reverse() # reverse the shape length coordinates (going backwards from were the construction of half the shape mesh ended)
        
        # for the Y dimension
        y_coord_reversed = y_coord.copy()
        y_coord_reversed = y_coord_reversed[0:-1]  # remove the last point
        y_coord_reversed = list(np.negative(y_coord_reversed))  # get the reflection of the half the shape mesh in the Y dimension
        y_coord_reversed.reverse()  # reverse the shape length coordinates (going backwards from were the construction of half the shape mesh ended)
        
        # Add the upper half and the lower half (reflection) of the shape mesh.
        x_coord = x_coord + x_coord_reversed
        y_coord = y_coord + y_coord_reversed
        
        # plot the shape mesh
        if show == True:
            plt.figure(figsize=(length,width))
            plt.plot(x_coord, y_coord)
    #        plt.legend()
            plt.xlabel('length (μm)',  size=10, fontweight='bold')
            plt.ylabel('width (μm)',  size=10, fontweight='bold')
            plt.show()
        
        # create a list with the shape mesh coordinates, compatible with the Polygon library
        shapely_points = []
        for c in range(len(x_coord)):
            shapely_points.append((x_coord[c], y_coord[c]))
            
        # create a shape polygon object
        polygon = Polygon(shapely_points)
        
        # a tuple with the x_coord and y_coord lists
        shape_mesh = (x_coord, y_coord)
        
        return shape_mesh, shapely_points, polygon
        
    
    
    def add_cell(self, cell_width, cell_length, verbose=False):
        """
        This function is used to add cell boundaries in the simulation.        
        A 2D_cell represents a cell region containing the cell boundaries.
        A cell has a width and a height.
        A hypothetical cell mesh is constructed using the specified width and the length.
        A bufferer zone is constructed on the inner side of the cell boundary, corresponding to the particle radius.
            The particle center cannot approach the cell boundary across this buffer zone (erosion) due to its size.
        
        Input:
            cell_width: positive float - the cell width in micrometers
            cell_length: positive float - the cell length in micrometers
        
        Returns:
            self.cell_mesh (x,y), self.cell_polygon: the cell mesh coordinates and the cell polygon are stored in the class
            self.cell_shapely_points: the POlygon dependent mesh coordinates ([x1,y1],[x2,y2],...) are also stored in the class
            self.added_cell variable is also updated to 1 from 0 to indicate that a cell was added and constrains the simulation accordingly
        """
        
        if self.added_cell == 1:
            print('warning: a cell has already been added. To proceed press Enter. To stop press ctl+C')
            input()
        
        if verbose == True:
            print('creating cell with length', cell_length, 'μm, and width', cell_width, 'μm...')
        
        # store the cell width and the cell length in the class
        self.cell_width = cell_width
        self.cell_length = cell_length
        
        # get the bacilus cell shape
        cell_mesh, shapely_points, cell_polygon = self.design_geometrical_shape(width=cell_width, length=cell_length)
        
        # store the cell polygon object in the class
        self.cell_polygon = cell_polygon
        
        # erode the cell polygon to simulate the exclusion of the particle center by the cell boundary
        # the amount of erosion is the same as the radius of the particle
        if self.particle_radius > 0: # if the particle radius is considered
            cell_polygon_eroded = cell_polygon.buffer(-self.particle_radius)
            self.cell_polygon_eroded = cell_polygon_eroded
        elif self.particle_radius == 0: # if the particle radius is negligible
            self.cell_polygon_eroded = self.cell_polygon
        
        # store the cell mesh coorinates in the class
        self.cell_shapely_points = shapely_points
        self.cell_mesh = cell_mesh
        
        # update the added_cell parameter since a cell was added.
        self.added_cell = 1
        
    
    
    def add_nucleoid(self, nucleoid_width, nucleoid_length, verbose=False):
        """
        This function is used to add a nucleoid when a cell object has already been specified.        
        A 2D_nucleoid represents a nucleoid region restricted within the cell boundaries.
        A nucleoid has a width and a height.
        A hypothetical nucleoid mesh is constructed using the specified width and the length.
        A bufferer zone is constructed on the outer side of the nucleoid boundary, corresponding to the particle radius.
            The particle center cannot approach the nucleoid boundary across this buffer zone (dilation) due to its size.
        
        Input:
            nucleoid_width: positive float - the nucleoid width in micrometers
            nucleoid_length: positive float - the nucleoid length in micrometers
        
        Returns:
            self.nucleoid_mesh (x,y), self.nucleoid_polygon: the nucleoid mesh coordinates and the nucleoid polygon are stored in the class
            self.nucleoid_shapely_points: the Polygon dependent mesh coordinates ([x1,y1],[x2,y2],...) are also stored in the class
            self.added_nucleoid variable is also updated to 1 from 0 to indicate that a nucleoid was added and constrains the simulation accordingly
        """
        
        class SimulationError(Exception):
            pass
        
        if self.added_cell == 0:
            raise SimulationError('No cell has been added to the simulation')
        elif self.added_cell == 1:
            if verbose == True:
                print('Creating nucleoid for cell with length',self.cell_length,'μm, and width',self.cell_width,'μm...')
                print('creating nucleoid with length', nucleoid_length, 'μm, and width', nucleoid_width, 'μm...')
        
        # store the nucleoid length and width in the class
        self.nucleoid_width = nucleoid_width
        self.nucleoid_length = nucleoid_length
    
    
        # get the nucleoid shape
        nucleoid_mesh, shapely_points, nucleoid_polygon = self.design_geometrical_shape(width=nucleoid_width, length=nucleoid_length)
        
        # store the cell polygon object in the class
        self.nucleoid_polygon = nucleoid_polygon
        
        # dilate the nucleoid polygon to simulate the exclusion of the particle center by the nucleoid boundary
        # the amount of dilation is the same as the radius of the particle
        if self.particle_radius > 0: # if the particle radius is considered
            nucleoid_polygon_dilated = nucleoid_polygon.buffer(self.particle_radius)
            self.nucleoid_polygon_dilated = nucleoid_polygon_dilated
        elif self.particle_radius == 0: # if the particle radius is negligible
            self.nucleoid_polygon_dilated = self.nucleoid_polygon
        
        # store the nucleoid mesh coorinates in the class
        self.nucleoid_shapely_points = shapely_points
        self.nucleoid_mesh = nucleoid_mesh
        
        # update the added_nucleoid parameter since a nucleoid was added to an existing cell.
        self.added_nucleoid = 1     
        
    
    
    def show_cell(self):
        
        if self.added_cell == 1:
            plt.figure(figsize=(self.cell_length,self.cell_width))
            plt.plot(self.cell_mesh[0], self.cell_mesh[1], color='blue')
            print('cell with length',self.cell_length,'μm, and width',self.cell_width,'μm...')
            plt.xlabel('length (μm)',  size=10, fontweight='bold')
            plt.ylabel('width (μm)',  size=10, fontweight='bold')
        else:
            print('no cell added')
        if self.added_nucleoid == 1:
            plt.plot(self.nucleoid_mesh[0], self.nucleoid_mesh[1], color='red')
            print('that contains a nucleoid with length',self.nucleoid_length,'μm, and width',self.nucleoid_width,'μm')
        else:
            print('no nucleoid added')
        plt.show()
        
    
    
    def get_cell_polygon(self):
        """
        Returns the cell polygon object
        """
        
        if self.added_cell == 1:
        # if the particle radius is 0 then the cell_polygon and the cell_polygon_eroded objects are the same
            return self.cell_polygon, self.cell_polygon_eroded
            
        elif self.added_cell == 0:
            print('No cell was aded.')
 
    
    
    def get_nucleoid_polygon(self):
        """
        Returns the nucleoid polygon object
        """
        
        if self.added_nucleoid == 1:
        # if the particle radius is 0 then the nucleoid_polygon and the nucleoid_polygon_dilated are the same
            return self.nucleoid_polygon, self.nucleoid_polygon_dilated
        
        elif self.added_nucleoid == 0:
            print('No nucleoid was added')
    
    
    
    def get_active_polygon(self):
        """
        Returns an active polygon within which a particle can move.
        
        """
        
        if self.added_cell == 1:  # this code works only if a cell was added
            
            if self.added_nucleoid == 0: # if a nucleoid was NOT added
                
                self.active_polygon = self.get_cell_polygon()[1] # the active polygon is the cell polygon (eroded or not)
                
                return self.active_polygon
                
            elif self.added_nucleoid == 1: # if a nucleoid was added
                
                cell_polygon = self.get_cell_polygon()[1] # get the cell polygon (eroded or not)
                nucleoid_polygon = self.get_nucleoid_polygon()[1] # get the nucleoid polygon (dilated or not)
                
                # get the symmetric difference between the nucleoid and the cell polygon and select the intersection between this symmetric difference and the cell polygon
                # This returns the coordinates that belong to the cell but not the nucleoid. These coordinates corespond to the active polygon
                active_polygon = (cell_polygon.symmetric_difference(nucleoid_polygon)).difference(nucleoid_polygon) 
                
                self.active_polygon = active_polygon
                
                return self.active_polygon
        
        elif self.added_cell == 0:
            print('No cell was added.')
        
            
    
    def is_position_in_cell(self, pos):
        """
        Checks if a selected position is in the cell

        pos: a position tuple (x,y).
        returns: True if pos is in the cell, False otherwise.
        """
        
        x, y = pos
        
        if self.added_cell == 1:
            
            return self.get_cell_polygon()[1].contains(Point(x,y))

        elif self.added_cell == 0:
            print('No cell was added.')
            
    
    
    def is_position_in_nucleoid(self, pos):
        """
        Checks if a selected position is in the nucleoid

        pos: a position tuple (x,y).
        returns: True if pos is in the nucleoid, False otherwise.
        """
        
        x, y = pos
        
        if self.added_nucleoid == 1:
            return self.get_nucleoid_polygon()[1].contains(Point(x,y))
        
        elif self.added_nucleoid == 0:
            print('No nucleoid added.')
            
    
    
    def is_position_in_active_polygon(self, pos):
        """
        Checks if a selected position is in the active polygon (polygon space that the particle can explore)

        pos: a position tuple (x,y).
        returns: True if pos is in the active polygon, False otherwise.
        """
        
        x, y = pos
        
        if self.added_cell == 1:
            return self.get_active_polygon().contains(Point(x,y))
        
        elif self.added_cell == 0:
            print('No cell was added.')

    
    
    def nucleoid_crossing(self, trajectory_trace):
        """
        Checks if a trajectory displacement crosses the nucleoid boundaries.
        
        Input:
            trajectory_trace: a list of two tuples corresponding to the coordinates of the trajectory trace during one displacement [(old_x, old_y),(new_x, new_y)]
        Returns:
            True or False if the particle displacement crosses the nucleoid or not.
        """
        
        nucleoid_polygon = self.get_nucleoid_polygon()[1]
        trajectory_line = LineString(trajectory_trace)
        
        return trajectory_line.intersects(nucleoid_polygon)
        
    
    
    def random_point_within_cell(self, num_points=1):
        """
        Select a number of points from within the cell that is not in the nucleoid (if a nucleoid is constructed).
        
        Input:
            num_points: positive integer - the number of points to be returned
            
        Returns:
            A list of the randomly selected points.
        """
        
        if self.added_cell == 1: # Check if a cell was adde
            active_polygon = self.get_active_polygon()
            
            min_x, min_y, max_x, max_y = active_polygon.bounds # get the minimum and maximum rectangular boundaries of the active polygon
            
            list_of_points = [] # initialize a list to store the randomly selected points (shapely Point objects)
            list_of_point_coordinates = []
    
            while len(list_of_points) < num_points: # select random points until the specified number of randomly selected points is reached (num_points)
                pnt = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y)) # select a random point from within the rectangular bounds of the active polygon
                if active_polygon.contains(pnt): # if the point belongs to the active polygon
                    list_of_points.append(pnt) # append the random point to the list (one point was added to the list)
                    list_of_point_coordinates.append(pnt.coords[0])
            self.random_points = list_of_points # store the randomly selected points in the class
            self.random_point_coordinates = list_of_point_coordinates
            
            return self.random_points, self.random_point_coordinates

        elif self.added_cell == 0:
            print('No cell was added. Cannot return a random cell position')
        
        
    
    def set_initial_position(self, pos):
        """
        Set the initial position of the particle in the simulation.
        This function can also be used to initialize a new round of simulation
        
        Input:
            pos: a tuple with the coordinates of the original position
        Returns:
            updates the displacement lists
        """
        
        if self.added_cell == 1: # if a cell was added
            point = Point(pos[0], pos[1])
            
            if self.get_active_polygon().contains(point):
                self.dx_list = []
                self.dy_list = []
                
                self.positions_x_list = []
                self.positions_y_list = []
                
    #            self.dx_list.append(pos[0])
    #            self.dy_list.append(pos[1])
                self.x = pos[0]
                self.y = pos[1]
                
                self.positions_x_list.append(pos[0])
                self.positions_y_list.append(pos[1])
            
            else:
                print('The specified point does not belong into the active polygon.')
        
        elif self.added_cell == 0: # if a cell was not added
            self.dx_list = []
            self.dy_list = []
            
            self.positions_x_list = []
            self.positions_y_list = []
            
            self.x = pos[0]
            self.y = pos[1]
            
            self.positions_x_list.append(pos[0])
            self.positions_y_list.append(pos[1])

    


    def get_displacement_list(self):
        """
        Get the list of x/y displacements until the current point in the simulation.
        """
        
        return self.dx_list, self.dy_list

    
    
    def get_position(self):
        """
        Get the current position of the particle
        """
        
        if len(self.positions_x_list) > 0: # if there were particle positions stored
        
            dx_list, dy_list = self.get_displacement_list()
            
            x = sum(dx_list) + self.positions_x_list[0]
            y = sum(dy_list) + self.positions_y_list[0]
            
            self.x = x
            self.y = y
            
            return self.x, self.y
        
        elif len(self.positions_x_list) == 0: # if there were not particle positions stored
            print('No particle positions were stored. Set the intial particle position')
    
    
    
    def get_trajectory_positions(self):
        """
        Returns two lists with all the particle positions in the x and y dimensions
        """
        
        return self.positions_x_list, self.positions_y_list
    


    def set_next_position(self, data):
        """
        Get the next position of the particle using random sampling from a distribution of displacements.
        
        Input:
            The diffusion coefficient in μm^2/sec (float) or a list/numpy-array with all the measured dispalcements
        Returns:
            Updates the displacements list including the new displacements randomly selected from a distribution of displacements.
                The displacements are only included in the displacement lists if they belong into the active polygon.
                This is a recursive algorithm. It runs the set_next_position function until a displacement that leads into a final position within the active polygon is selected.
                Also those displacements that cross the nucleoid are excluded.
        """
        
        if type(data) == float or type(data) == int: # if the data is the diffusion coefficient (signle float value)
            
            diffusion_coefficient = data
            
            std_normal = np.sqrt(2*diffusion_coefficient*self.interval) # The function which gives the std of the normal distribution from the Diffusion coefficient
        
            dx_random = np.random.normal(0, std_normal, 1)[0]
            dy_random = np.random.normal(0, std_normal, 1)[0]
        
        elif type(data) == list or type(data) == np.ndarray:
            if len(data) > 1: # if the data is a list or numpy array with the experimentally measured displacements
                distribution_of_displacements = data
                
                dx_random = random.sample(list(distribution_of_displacements),1)[0]
                dy_random = random.sample(list(distribution_of_displacements),1)[0]
        
        
        
        old_x, old_y = self.get_position() # get the current particle position (before displacement)
        dx_list_old, dy_list_old = self.get_displacement_list() # get the current lists of all displacements
        x_list_old, y_list_old = self.get_trajectory_positions() # get all the positions of the trajectory so far
        
        dx_list_new = dx_list_old.copy() # create copies of the displacement lists
        dy_list_new = dy_list_old.copy()
        dx_list_new.append(dx_random) # update the new list of displacements
        dy_list_new.append(dy_random)
        
        new_x = sum(dx_list_new) + x_list_old[0] # get the new particle position
        new_y = sum(dy_list_new) + y_list_old[0]
        
        if self.added_cell == 0:  # if no cell has been added (simulation in unconstrained space)
            self.dx_list =  dx_list_new
            self.dy_list =  dy_list_new
            self.positions_x_list.append(new_x)
            self.positions_y_list.append(new_y)
        
        elif self.added_cell == 1: # if a cell has been added (cell confinement included)
            
            if self.is_position_in_active_polygon(pos=(new_x, new_y)) == True: # if the new position belongs in active polygon
                
                if self.added_nucleoid == 0: # if a nucleoid was NOT added
                    # since the particle belongs to the active polygon then the new particle position is returned
                    self.dx_list =  dx_list_new
                    self.dy_list =  dy_list_new
                    self.positions_x_list.append(new_x)
                    self.positions_y_list.append(new_y)
                    
                elif self.added_nucleoid == 1: # if a nucleoid was added things are more complicated
                    # even if the new particle position falls within the active polygon it may still cross the nucleoid
                    if self.nucleoid_crossing(trajectory_trace=[(old_x, old_y),(new_x, new_y)]) == False: # if the particle trajectory does not cross the nucleoid
                        self.dx_list =  dx_list_new
                        self.dy_list =  dy_list_new
                        self.positions_x_list.append(new_x)
                        self.positions_y_list.append(new_y)
                    elif self.nucleoid_crossing(trajectory_trace=[(old_x, old_y),(new_x, new_y)]) == True:
                        self.set_next_position(data)
            
            elif self.is_position_in_active_polygon(pos=(new_x, new_y)) == False: # if the new position does not belong in active polygon
                self.set_next_position(data)
            
    
    
    def show_particle(self, show, save_path = None):
        """
        Shows the particle trajectory.
        If a cell is added, the cell boundaries are also displayed.
        If a nucleoid is added, the nucleoid boundaries are also displayed.
        
        The first position of the trajectory is marked as green and the last one as red.
        
        The dashed lines indicate the dilated nucleoid boundaries or the eroded cell boundaries.
        The amount of dilation or erosion is equal to the particle radius, to illustrate the exclusion of the particle center from the cell boundaries or the nucleoid mesh.
        """
        
        if self.added_cell == 1:
            plt.figure(figsize=(self.cell_length*5, self.cell_width*5))
            cell_polygon = self.get_cell_polygon()
            x_cell, y_cell = cell_polygon[0].exterior.coords.xy
            x_cell_eroded, y_cell_eroded = cell_polygon[1].exterior.coords.xy
            
            plt.plot(x_cell,y_cell, color='pink', linewidth=3)
            plt.plot(x_cell_eroded,y_cell_eroded, '--', color='red', linewidth=0.5)
            
            if self.added_nucleoid == 1:
                nucleoid_polygon = self.get_nucleoid_polygon()
                x_nuc, y_nuc = nucleoid_polygon[0].exterior.coords.xy
                x_nuc_dilated, y_nuc_dilated = nucleoid_polygon[1].exterior.coords.xy
                
                plt.plot(x_nuc,y_nuc, color='green', linewidth=3)
                plt.plot(x_nuc_dilated,y_nuc_dilated, '--', color='turquoise', linewidth=0.5)
            
        
        plt.plot(self.positions_x_list,self.positions_y_list, color='silver', linewidth=2)
        plt.plot(self.positions_x_list[0],  self.positions_y_list[0], 'o', color='green', markerfacecolor='none', markersize=10)
        plt.plot(self.positions_x_list[-1], self.positions_y_list[-1], 'o', color='red', markerfacecolor='none', markersize=10)
        
        
        plt.xlabel('length (μm)', size=10, fontweight='bold')
        plt.ylabel('width (μm)', size=10, fontweight='bold')
        
        plt.axvline(0, linewidth=0.1, color='black')
        plt.axhline(0, linewidth=0.1, color='black')
        
        plt.text(0,0, s=('time: '+format(((len(self.positions_x_list)-1)*self.interval), '.2f')+' sec'), color='moccasin', fontsize=30, fontweight='bold', alpha=0.8)
        
        plt.box(False)
        
        if save_path != None:
            plt.savefig(save_path)

        if show == True:
            plt.show()
        else:
            plt.close()



def run_simulation(number_of_particles, cell_length, cell_width, nucleoid_width, 
                   nucleoid_length, time_steps, diffusion_coefficient_array, time_lag, 
                   particle_radius, random_seeds, nucleoid_confinement, show_interval, save_path):
    """
    Runs the simulation
    """
    
#    number_of_particles=1
#    length = 3 
#    width = 0.8
#    time_steps = 200
#    diffusion_coefficient = 0.1
#    time_lag = 0.5
#    scale = 5
#    particle_radius = 0.02
#    nucleoid_width=0.5
#    nucleoid_length=1.7
    
    np.random.seed(random_seeds)
    random.seed(random_seeds)
    
    particle_displacements = {}
    particle_positions = {}
    particle_msd = {}
    
    mean_dc = np.mean(np.log10(diffusion_coefficient_array))
    std_dc = np.std(np.log10(diffusion_coefficient_array))

    print('generating trajectories from a normal distribution of effective')
    print('diffusion coefficients with a mean of', 10**mean_dc, 'μm2/sec')

    sampled_diffusion_coefficients = 10**np.random.normal(mean_dc, std_dc, number_of_particles)
    
    for particle in range(number_of_particles):
        
        sim = bound_particle_simulation(particle_radius=particle_radius, interval=time_lag)
        
        # if a cell is added use this
        sim.add_cell(cell_width, cell_length)
        if nucleoid_confinement == True:
            sim.add_nucleoid(nucleoid_width, nucleoid_length)
        if particle == 0:
            sim.show_cell()
        point = sim.random_point_within_cell()[1][0]
        sim.set_initial_position(point)
        
        #if not cell us added use this
        # sim.set_initial_position((0,0))
        
#        sim.show_particle(save_path=save_path)

        for iteration in range(time_steps-1):
            
#             print(iteration+1, 'out of', time_steps)
#             print(sampled_diffusion_coefficients[particle])
            sim.set_next_position(float(sampled_diffusion_coefficients[particle]))
            
        if particle%show_interval==0:
            if save_path != None:
                sim.show_particle(False, save_path=save_path+'/particle'+str(particle)+'.jpeg')
            else:
                sim.show_particle(True, None)
        
        particle_displacements[particle] = sim.get_displacement_list()
        particle_positions[particle] = sim.get_trajectory_positions()
#         particle_msd[particle] = sim.msd_estimation()[0]


    return particle_displacements, particle_positions