import numpy as np
from netCDF4 import Dataset

nx = 250
ny = 1000
n_layers = 40 
n_cells = nx * ny
n_edges = ( n_cells * 3 ) + (nx * 2)
rho = 1026.0
ncfile_input = Dataset('initial_output.nc', 'r')
ncfile_output = Dataset('velocity_init.nc', 'r+')

variables = ( 
      {"name": "pressure", "var": "pressure", "y_val":[]},
      {"name":"coriolis", "var": "fCell", "y_val":[]},
      {"name": "y", "var": "yCell", "y_val":[]})

velocity_full_file = ncfile_output.variables["normalVelocity"]


for variable in variables:
    var_full = [] 
    # Extract variable fro output file
    var_full_file = ncfile_input.variables[variable["var"]]
    # Create a nested list for pressure by layers
    if variable["name"] == "pressure":
        for i in range(0, n_layers):
            var_full_reshape = np.reshape(var_full_file[:, :,i], [n_cells])
      	    var_full_layer = var_full_reshape.tolist()
            var_full.append(var_full_layer)
    # Create a list for the other single layer variables
    else:
        var_full_reshape = np.reshape(var_full_file[:], [n_cells])
        var_full_layer = var_full_reshape.tolist()
        var_full.append(var_full_layer)

    # Calculate variables at midpoints and create a list in terms of y for each layer
    for i in range(0, len(var_full)):
        # Pull out one value from each row
        var_center = var_full[i][0::nx]
        # Initialize midpoint array with zeros
        var_mid = [0] * (ny + 1)
        # Put center and midpoint arrays into main array
        var_y_layer = var_center + var_mid
        var_y_layer[1::2] = var_center
        var_y_layer[0::2] = var_mid
        # Find midpoint values by averaging neighboring values
        for i in range(2, len(var_y_layer) - 1, 2):
            var_y_layer[i] = ( var_y_layer[i-1] + var_y_layer[i+1] ) / 2  
        # Find yCell values at the top and bottom edges
        if variable["name"] == "y":
            var_y_layer[0] = var_y_layer[1] - (var_y_layer[2] - var_y_layer[1])
            var_y_layer[-1] = var_y_layer[-2] + (var_y_layer[2] - var_y_layer[1])
        # Add each layer to the nested list for the full value
        variable["y_val"].append(var_y_layer)

velocity_full = []
# Calculate velocity at y locations and map to edges for each layer
for j in range(0, n_layers):
    # Initialize velocity arrays
    velocity_y_layer = [0] * ((ny * 2) + 1)
    velocity_edges = [0] * n_edges
    velocity_edges_center = [0] * n_cells
    velocity_edges_left = [0] * (n_cells + nx)
    velocity_edges_right = [0] * (n_cells + nx)
    # Calculate velocity
    for i in range(2, len(velocity_y_layer)-2):
        velocity_y_layer[i] = ( (-1)/(variables[1]["y_val"][0][i] * rho) ) * \
        ( (variables[0]["y_val"][j][i+1] - variables[0]["y_val"][j][i-1]) / \
        (variables[2]["y_val"][0][i+1] - variables[2]["y_val"][0][i-1]) )
    # Make velocities at top and bottom cells 0.5 of neighboring cells
    velocity_y_layer[1] = 0.5 * velocity_y_layer[2]
    velocity_y_layer[-2] = 0.5 * velocity_y_layer[-3]
    # Map values from y array to edge arrays
    # Calculate velocity in y direction at midpoint edges
    velocity_edges_center = np.repeat(velocity_y_layer[1::2], nx)

    velocity_edges_left = np.repeat(velocity_y_layer[0::2], nx)
    for i in range(0, len(velocity_edges_left)):
        velocity_edges_left[i] = 0.5 * velocity_edges_left[i]

    velocity_edges_right = np.repeat(velocity_y_layer[0::2], nx)
    for i in range(0, len(velocity_edges_right)):
        velocity_edges_right[i] = -0.5 * velocity_edges_right[i]
    # Create full edge array
    top_length = nx * 2
    velocity_edges[0:n_edges-top_length:3] = velocity_edges_center
    velocity_edges[1:n_edges-top_length:3] = velocity_edges_left[0:-nx]
    velocity_edges[2:n_edges-top_length:3] = velocity_edges_right[0:-nx]

    velocity_full_file[0, :, j] = velocity_edges    
    # Add layer to full depth edge array
    #velocity_full.append(velocity_edges)
ncfile_input.close()
ncfile_output.close()    

    


