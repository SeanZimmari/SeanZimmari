#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np

import sympy as sp

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.interpolate import interp1d

from sklearn.metrics import mean_squared_error

def mesh(Nelem, L):
    
    #This function will split the bar into (2*Nelem) + 1 nodes
    
    z_domain = np.linspace(0, L, (2 * Nelem) + 1);
    
    return z_domain
    

def Assembly_shape_functions(Nelem, Length_element):
    
    #This function will assign to each node the specific shape function
    #However, I will consider for each element equal shape functions, but I will have to consider
    #in the relationship the length of each element Lenght_element
    
    #Initialization of a symbolic z variable
    
    z_symbol = sp.symbols('z');
    
    #Initialization of the shape_functions ditionary
    
    shape_functions = {};
    
    #Initialization of the keys -----> 'Node 1', 'Node 2', ....
    
    keys = [];
    
    #Initialization of the values -----> shape_functions
    
    values = [];
    
    #Initialization of the variable which will take into account the number of nodes
    
    element = 0;
    
    #Initialization of a list that will contain the positions of the nodes for each element
    
    nodes_element = [];
    
    #Initialization of a variable which will take into account the number of iterations made for the nodes
    
    iterations = 0
    
    #Initialization of a variable which will take into account the different lengths
    
    node = 0;
    
    #We will have to iterate for every element
    
    for element in range(Nelem):
        
        z_i = 0;
        
        z_j = Length_element/2;
        
        z_k = Length_element
        
        Ni = (2/Length_element**2) * (z_symbol - z_j)* (z_symbol - z_k);

        node += 1;

        keys.append(f'Node {node} Element {element + 1}')

        values.append(Ni);

        Nj = -(4/Length_element**2) * (z_symbol - z_i)* (z_symbol - z_k);

        node += 1;

        keys.append(f'Node {node} Element {element + 1}')

        values.append(Nj);

        Nk = (2/Length_element**2) * (z_symbol - z_i)* (z_symbol - z_j);

        node += 1;

        keys.append(f'Node {node} Element {element + 1}')

        values.append(Nk);
        
        node = 0;
    
    #Assembly the shape_functions dictionary
    
    for key, value in zip(keys,values):
        
        shape_functions[key] = value;
        
    return shape_functions
    
def compute_derivatives_of_shape_functions(shape_functions):
    
    #This function will derive the shape functions
    
    #Definition of z symbol
    
    z_symbol = sp.symbols('z');
    
    #Definition of the array that will contain the derivatives of the shape functions
    
    derivatives_shape_functions = {};
    
    #Definition of keys and values
    
    keys = [];
    
    values = [];
    
    #Define a variable which will take into account the mapping in the shape functions
    
    shape = 0;
    
    element = 0;
    
    #Derivation of shape functions
    
    for function in range(int(len(shape_functions))):
        
        if shape < 3:
        
            derivative = (sp.diff(shape_functions[f'Node {shape + 1} Element {element + 1}'], z_symbol));
            
            values.append(derivative);
            
            keys.append(f'Node {shape + 1} Element {element + 1}');

            shape += 1;
        
        else:
            
            shape = 0;
            
            element += 1;
            
            derivative = (sp.diff(shape_functions[f'Node {shape + 1} Element {element + 1}'], z_symbol));
            
            values.append(derivative);
            
            keys.append(f'Node {shape + 1} Element {element + 1}');

            shape += 1;
            
    #Assembly the derivatives_shape_functions dictionary
    
    for key, value in zip(keys,values):
        
        derivatives_shape_functions[key] = value;       
    
    return derivatives_shape_functions

                
def local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A):

    #Initialize the local stiffness matrixes. This variable will contain all the local stiffness matrixes
    
    local_stiffness_matrixes = [];

    #Initialize a symbolical variable

    z_symbol = sp.symbols('z')

    #Initialize a list which will contain the derivatives of the shape functions

    #Get the derivatives of the shape functions

    for element in range(Nelem):
        
        #Initialize a list which will contain the derivatives of the shape functions
        
        dN = [];
        
        #Get the shape functions for that element: the first element will be equal to the shape function at the node 2*i
        #first iteration it will get the first shape function, second iteration it will get the third shape_function and so on

        for i in range(3):

            dN.append(derivatives_shape_functions[f'Node {i + 1} Element {element + 1}'])

        #Hence, once we have initialized the element stiffness, K[i, j] = EA/L * integral(dNi/dx * dNj/dx)^2. I will set
        #the integral between 0 to the length of the element z

        element_stiffness = np.zeros([3, 3]);

        for i in range(3):

            for j in range(3):

                integrand = dN[i] * dN[j];

                element_stiffness[i, j] = (E * A) * sp.integrate(integrand, (z_symbol, 0, Length_element));

        local_stiffness_matrixes.append(element_stiffness);      

    return np.array(local_stiffness_matrixes)

def Assembly_connectivity_matrix(Nelem):
    
    #This function will assembly the connectivity matrix
    #Every element has 3 nodes so the connectivity matrix will have 3 columns and Nelem rows
    
    Connectivity_matrix = np.zeros([Nelem, 3]);
    
    for i in range(Nelem):
        
        #At the first iteration make the three elements
        
        if i == 0:
        
            for j in range(3):

                    Connectivity_matrix[i,j] = j + 1;

                    j += 1;
                
        else:
            
            #At higher iterations, the first element will be equal to the one at the previous row and third column
            
            Connectivity_matrix[i, 0] = Connectivity_matrix[i - 1, 2];
            
            for j in range(2):

                Connectivity_matrix[i,j + 1] = Connectivity_matrix[i - 1, 2] + j + 1;

                j += 1;
        
        i += 1;
        
    return Connectivity_matrix

def Assembly_global_stiffness_matrix(Nelem, Connectivity_matrix, Total_Element_stiffness_matrix):
    
    K_g = np.zeros([(2*Nelem) + 1, (2*Nelem) + 1]);
    
    for index in range(Nelem):
    
        stiffness_matrix = np.zeros([3, 3]);

        element_connectivity_matrix = np.zeros(3);
        
        stiffness_matrix = Total_Element_stiffness_matrix[index];
        
        element_connectivity_matrix = Connectivity_matrix[index, :];
        
        for i in range(3):
            
            for j in range(3):
                
                K_g[int(element_connectivity_matrix[i] - 1), int(element_connectivity_matrix[j] - 1)] += stiffness_matrix[i][j];
                
                j += 1;
                
            i += 1;
        
        index += 1
                
    return K_g

def Assembly_force_vector(Nelem, Length_element, a, b):
    
    #Define the force vector: it will have the same size as the number of nodes
    
    F_global = np.zeros([(2 * Nelem) + 1]);
    
    #Iterate between the different elements
    
    for element in range(Nelem):
    
        #Define the new a (a is gonna vary between the different elements):

        a_new = a + b*Length_element*element;
        
        #Define the positions of the nodes: with this notations, at the beginning of each iteration we will come back to
        #the node of the previous iteration
        
        position_node = [2 * element, (2 * element) + 1, (2 * element) + 2];
        
        #Definition of the forces(Same relationship of the analytical solution, we have only to substitute the different values of z)
        #NOTE: Obviously we have to consider the new a and the length of the element
        
        Local_forces = np.array([(a_new*Length_element/6), (b*Length_element**2/3 + 2*a_new*Length_element/3), (b*Length_element**2/6 + a_new*Length_element/6)]);
        
        #Place the forces of each element into a global vector
        
        for i in range(3):
            
            F_global[position_node[i]] += Local_forces[i];
            
            i += 1;
            
        element += 1;
        
    return F_global

def apply_boundary_conditions(K_G, BC, F_global): # Setting buondary condition for zero and non-zero displacements

        K_mod = K_G.copy() #copies for manipulation
        F_mod = F_global.copy()

        nodes_with_disp = BC[0] - 1  # Adjust indexing
        displacements = BC[1]

      
        sorted_indices = np.argsort(nodes_with_disp)[::-1]   # Sort nodes in descending order for index shifting
        nodes_with_disp = nodes_with_disp[sorted_indices]
        displacements = displacements[sorted_indices]

        colrem = np.zeros(len(F_mod)) #Column reminder of k equivalent stiffness elements

        for i, node in enumerate(nodes_with_disp): # Deleting raws and columns according to BC, plus  subtracting the k stiffness elements reminder column to force vector 
            disp = displacements[i]
            if disp == 0:# Nodes with zero displacement
                K_mod = np.delete(K_mod, node, axis=0)
                K_mod = np.delete(K_mod, node, axis=1)
                F_mod = np.delete(F_mod, node)
                colrem = np.delete(colrem, node)
            else: #Nodes with non-zero displacement
                if node < K_mod.shape[1]:#  Ensuring that the node still corresponds to a valid column in updated K_G
                    colrem[:len(K_mod)] -= K_mod[:, node] * disp
                K_mod = np.delete(K_mod, node, axis=0)
                K_mod = np.delete(K_mod, node, axis=1)
                F_mod = np.delete(F_mod, node)
                colrem = np.delete(colrem, node)

        result = F_mod + colrem[:len(F_mod)]
        return K_mod, result
    
def Computation_of_displacements(New_global_stiffness, New_force_vector, Nelem, BC):
    
    #This function will compute the displacements after modifications are applied to the global stiffness and to the
    #Force vector
    
    displacements = np.zeros([(2 * Nelem) + 1]);
    
    #Debugging
    
    if New_global_stiffness.size > 0:
        
        #Solve the system
        
        Finaldisplacements = np.linalg.solve(New_global_stiffness, New_force_vector)
        
    else:
        
        Finaldisplacements = np.array([])  # If BC applied to all the nodes then no displacements to find

   
    for i, node in enumerate(BC[0]):# Add to the total displacements vector the displacements of the nodes with BC by re-indexing
        
        displacements[node - 1] = BC[1][i]

    
    free_nodes = [i for i in range(2*Nelem + 1) if i not in (BC[0] - 1)] # Free nodes indexing from BC, then add to the total displacements vector the found free displacements
    
    for i, node in enumerate(free_nodes):
        
        displacements[node] = Finaldisplacements[i]
    
    return displacements

def scaling_displacements(z_domain, Nelem, shape_functions, displacements, Length_element):
    
    #This function will scale the displacements along the length of the bar
    
    z_symbol = sp.symbols('z')  # Ensure 'z' is defined as a symbolic variable
    
    plotting_points_per_element = 50  # Number of refined points for each element
    
    # Initialize lists for the scaled displacements and corresponding points in the domain
    
    scaled_displacements = []
    
    scaled_points = []
    
    # Iterate over each element to compute the refined displacements
    
    for elem in range(Nelem):
        # Get nodal displacements for the current element
        
        disp_node1 = displacements[2 * elem]          # Displacement at the first node of the element
        
        disp_node2 = displacements[2 * elem + 1]      # Displacement at the middle node
        
        disp_node3 = displacements[2 * elem + 2]      # Displacement at the third node
        
        # Nodal displacements array
        disp = [disp_node1, disp_node2, disp_node3]
        
        # Define local plotting points for the current element
        local_plotting_points = np.linspace(0, Length_element, plotting_points_per_element)
        
        # Compute the displacements at the refined plotting points using shape functions
        for point in local_plotting_points:
            numerical_function = 0
            
            # Multiply nodal displacements by shape functions evaluated at the refined point
            for j in range(3):
                numerical_function += disp[j] * shape_functions[f'Node {j + 1} Element {elem + 1}'].subs(z_symbol, point).evalf()
            
            # Append the refined displacement and corresponding point in the domain
            scaled_displacements.append(float(numerical_function))
            
            scaled_points.append(point + (elem * Length_element))
    
    return scaled_displacements, scaled_points

def Compute_stress_strain(scaled_displacements, Nelem, E):
    
    #This function will compute the strains and stresses
    
    strain = np.zeros([Nelem]);
    
    stress = np.zeros([Nelem]);
    
    #Computation of the strains
    
    for i in range(Nelem):
        
        strain[i] = (scaled_displacements[(2 * i) + 2] - scaled_displacements[2 * i])/Length_element;
        
        i += 1;
    
    #Computation stresses
    
    for j in range(Nelem):
        
        stress[j] = E * strain[j];
        
        j += 1;
        
    return strain, stress

def scale_stresses(Nelem, shape_functions, displacements, Length_element, E):
    
    # Function to compute and scale stresses along the length of the bar
    
    z_symbol = sp.symbols('z')  # Symbolic variable for position along the element
    
    element_stress = []   # List to store stress values for each element
    
    element_centers = []  # List to store the center positions of each element
    
    # Iterate through each element
    for elem in range(Nelem):
        
        # Get the nodal displacements for the current element (3 nodes per element)
        
        disp_node1 = displacements[2 * elem]      # Displacement at node 1
        
        disp_node2 = displacements[2 * elem + 1]  # Displacement at node 2
        
        disp_node3 = displacements[2 * elem + 2]  # Displacement at node 3
        
        # Combine displacements into a list
        
        displacements_elem = [disp_node1, disp_node2, disp_node3]
        
        # Compute the strain at the midpoints of the element
        for z in np.linspace(0, Length_element, 10):  # 10 points along the element
            
            strain_at_z = 0
            
            # Compute the shape function derivatives at each z
            
            for j in range(3):
                
                shape_func_derivative = sp.diff(shape_functions[f'Node {j + 1} Element {elem + 1}'], z_symbol)
                
                strain_at_z += displacements_elem[j] * shape_func_derivative.subs(z_symbol, z).evalf()
            
            # Compute stress using Hooke's law: stress = E * strain
            
            stress_at_z = E * float(strain_at_z)
            
            # Store the computed stress and position
            
            element_stress.append(stress_at_z)
            
            element_centers.append(z + (elem * Length_element))  # Store position along the length of the bar
    
    return element_stress, element_centers

# Function to calculate MSE error

def calculate_L2_error(fem_values, analytical_values):
    
    #This function will return the mean square error in order to evaluate the convergence rate
    
    return np.mean((np.array(analytical_values) - np.array(fem_values))**2)

def compute_mse_disp(L, numerical_analytical_values_zdomain, scaled_FEM_displacements, plotting_points, Nelem):
    
    delta_z = np.linspace(0, L, 400)  # Create a fine grid for evaluating analytical solutions.
    
    analytical_values_deltaz = np.array([numerical_analytical_values_zdomain(z) for z in delta_z])  # Evaluate analytical solutions.
    
    # Interpolate FEM displacements to the same fine grid.
    
    fem_interpolator = interp1d(plotting_points, scaled_FEM_displacements, kind='linear', fill_value="extrapolate")
    
    fem_values_deltaz = fem_interpolator(delta_z)
    
    # Initialize MSE for each element.
    
    mse_disp = np.zeros(Nelem)
    
    length_per_element = L / Nelem  # Calculate length per element.
    
    # Calculate MSE for each element.
    
    for elem in range(Nelem):
        
        z_start = elem * length_per_element
        
        z_end = (elem + 1) * length_per_element
        
        # Get points in delta_z that correspond to this element's domain.
        
        mask = (delta_z >= z_start) & (delta_z <= z_end)
        
        # Select analytical and FEM values for the current element.
        
        analytical_elem = analytical_values_deltaz[mask]
        
        fem_elem = fem_values_deltaz[mask]
        
        # Check for empty arrays and matching lengths
        
        if analytical_elem.size > 0 and fem_elem.size > 0 and len(analytical_elem) == len(fem_elem):
            
            mse_disp[elem] = mean_squared_error(analytical_elem, fem_elem)
            
        else:
            
            mse_disp[elem] = np.nan  # Assign NaN if there's an issue
    
    return mse_disp

def compute_mse_stress(L, Nelem, numerical_analytical_function, scaled_FEM_stresses):
    
    #This function will evaluate the mse for the stresses
    
    delta_z = np.linspace(0, L, 10*Nelem); #Discretizing points
    
    analytical_values_deltaz = np.zeros([len(delta_z)]);
    
    for point in delta_z:
        
        analytical_stress_deltaz = numerical_analytical_function(point);
    
    return calculate_L2_error(scaled_FEM_stresses, analytical_stress_deltaz)

def plot_L2_error_vs_elements(Nelem_list, mse_results):
    """
    Plots the L2 error vs number of elements for FEM and analytical solutions.

    Parameters:
    - Nelem_list: List of number of elements.
    - mse_results: List of MSE arrays (one for each Nelem).
    """
    
    # Extract mean square errors for plotting
    
    mean_square_errors = [np.mean(mse) for mse in mse_results]  # Average MSE over elements

    # Reverse the order if needed to ensure higher MSE corresponds to lower Nelem
    
    Nelem_list = np.array(Nelem_list)
    
    mean_square_errors = np.array(mean_square_errors)

    # Ensure we are plotting in decreasing order
    
    sorted_indices = np.argsort(Nelem_list)
    
    sorted_Nelem_list = Nelem_list[sorted_indices]
    
    sorted_errors = mean_square_errors[sorted_indices]

    # Plot the L2 error vs number of elements
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_Nelem_list, sorted_errors, marker='o', linestyle='-', color='b', label='MSE Error')
    plt.xlabel('Number of Elements (Nelem)')
    plt.ylabel('Mean Square Error Displacements')
    plt.title('Mean Square Error vs Number of Elements')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()
    
def compute_mse_stress(L, numerical_analytical_function, scaled_FEM_stresses, plotting_points, Nelem):
    """
    This function evaluates the mean square error for the stresses element-wise.
    
    Parameters:
    - L: Length of the bar.
    - numerical_analytical_function: Function that returns the analytical stress at any z-value.
    - scaled_FEM_stresses: Array of FEM stresses at plotting points.
    - plotting_points: Points in the bar where FEM stresses are evaluated.
    - Nelem: Number of elements used in the FEM.
    
    Returns:
    - mse_stress: Array of MSE values, one per element.
    """
    
    # Create a fine grid for evaluating analytical solutions.
    delta_z = np.linspace(0, L, 400)
    
    # Evaluate analytical solutions.
    analytical_values_stress = np.array([numerical_analytical_function(z) for z in delta_z])
    
    # Interpolate FEM stresses to the same fine grid.
    fem_interpolator = interp1d(plotting_points, scaled_FEM_stresses, kind='linear', fill_value="extrapolate")
    fem_values_stress = fem_interpolator(delta_z)
    
    # Initialize MSE for each element.
    mse_stress = np.zeros(Nelem)
    
    # Calculate length per element.
    length_per_element = L / Nelem
    
    # Calculate MSE for each element.
    for elem in range(Nelem):
        
        z_start = elem * length_per_element
        
        z_end = (elem + 1) * length_per_element
        
        # Get points in delta_z that correspond to this element's domain.
        
        mask = (delta_z >= z_start) & (delta_z <= z_end)
        
        # Select analytical and FEM values for the current element.
        
        analytical_elem = analytical_values_stress[mask]
        
        fem_elem = fem_values_stress[mask]
        
        # Calculate MSE for the current element.
        
        mse_stress[elem] = mean_squared_error(analytical_elem, fem_elem)
    
    return mse_stress

def plot_stress_error_vs_elements(Nelem_list, mse_results_stress):
    """
    Plots the L2 error for stresses vs number of elements for FEM and analytical solutions.

    Parameters:
    - Nelem_list: List of number of elements.
    - mse_results_stress: List of MSE arrays for stresses (one for each Nelem).
    """
    
    # Extract mean square errors for plotting
    mean_square_errors_stress = [np.mean(mse) for mse in mse_results_stress]

    # Reverse the order if needed to ensure higher MSE corresponds to lower Nelem
    Nelem_list = np.array(Nelem_list)
    
    mean_square_errors_stress = np.array(mean_square_errors_stress)

    # Ensure we are plotting in decreasing order
    sorted_indices = np.argsort(Nelem_list)
    
    sorted_Nelem_list = Nelem_list[sorted_indices]
    
    sorted_errors_stress = mean_square_errors_stress[sorted_indices]

    # Plot the L2 error vs number of elements for stresses
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_Nelem_list, sorted_errors_stress, marker='o', linestyle='-', color='r', label='Stress MSE Error')
    plt.xlabel('Number of Elements (Nelem)')
    plt.ylabel('Mean Square Error Stresses')
    plt.title('Mean Square Error vs Number of Elements')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()
    
#INPUT

A = 120 #mm^2
E = 70e3 #MPa
L = 500 #mm
a = 13 #N/mm
b = 0.13 #N/mm^2

#Define a library in order to plot the curves varying the number of element

Plotting_element = [1, 2, 3, 4, 5];

#Boundary condition

BC = np.array([[1],[0]]);

#Symbolical z for evaluation of analytical relationship

z_symbol = sp.symbols('z');

# Initialize lists to store displacement and stress results for each refinement

fem_displacement_results = []

fem_displacement_domain_results = []

fem_stress_results = []

fem_stress_domain_results = []

#Calculating the analytical displacements in z_domain

function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a*z_symbol**2)/2 - (a*L + (b*L**2)/2)*z_symbol))

numerical_function_analytical_displacement = sp.lambdify(z_symbol, function_analytical_displacement)

# Stress: σ(z) 

function_analytical_stress_nodes = -((1 / A) * ((b*z_symbol**2)/2 + (a * z_symbol) - (a * L + (b*L**2)/2)));

numerical_function_analytical_stress = sp.lambdify(z_symbol, function_analytical_stress_nodes)

mse_results_disp = [];

mse_results_stress = []

# Loop through different element sizes

for Nelem in Plotting_element:
    
    # Compute the FEM solution for the given number of elements
    
    Length_element = L / Nelem

    z_domain = mesh(Nelem, L)
    
    #Calculating the analytical displacements in z_domain

    function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a*z_symbol**2)/2 - (a*L + (b*L**2)/2)*z_symbol))

    numerical_function_analytical_displacement = sp.lambdify(z_symbol, function_analytical_displacement)
    
    analytical_displacements_zdomain = numerical_function_analytical_displacement(z_domain);
    
    #Computing shape functions
    
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    
    #Compute derivatives of shape functions
    
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    
    #Assembly element stiffness matrix
    
    stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    
    #Assembly connectivity matrix
    
    connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    
    #Assembly global stiffness matrix
    
    global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)
    
    #Compute force vector
    
    force_vector = Assembly_force_vector(Nelem, Length_element, a, b)
    
    #Apply boundary conditions
    
    new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

    # Compute FEM displacements
    
    FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)
    
    #Computing the L2_error

    # Scale displacements for plotting
    
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)
    
    # Compute FEM stresses
    
    FEM_strain, FEM_stresses = Compute_stress_strain(FEM_displacements, Nelem, E)
    
    #Scale stresses for plotting
    
    scaled_stresses, scaled_stress_domain = scale_stresses(Nelem, shape_functions, FEM_displacements, Length_element, E)
    
    # Compute MSE for displacements and stresses
    
    mse_disp = compute_mse_disp(L, numerical_function_analytical_displacement, scaled_displacements, scaled_points, Nelem);
    
    mse_results_disp.append(mse_disp)
    
    mse_stress = compute_mse_stress(L, numerical_function_analytical_stress, scaled_stresses, scaled_stress_domain, Nelem)
    
    mse_results_stress.append(mse_stress)
    
    # Store FEM displacement and stress results for plotting later
    
    fem_displacement_results.append(scaled_displacements)
    
    fem_displacement_domain_results.append(scaled_points)
    
    fem_stress_results.append(scaled_stresses)
    
    fem_stress_domain_results.append(scaled_stress_domain)

#Transform the mse in arrays for plotting later

mse_disp_results = np.array(mse_results_disp);

mse_stress_results = np.array(mse_results_stress);

# Generate analytical values over the domain

analytical_domain = np.linspace(0, L, 500)  # Fine grid for smooth plotting

analytical_displacements = numerical_function_analytical_displacement(analytical_domain)

analytical_stresses = numerical_function_analytical_stress(analytical_domain)

# Plot FEM and analytical solutions for displacements and stresses

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot Displacements

ax1.set_title('FEM Displacements vs Analytical Solution')

ax1.set_xlabel('Position along the bar (z) [mm]')

ax1.set_ylabel('Displacement [mm]')

ax1.grid(True)

# Plot FEM displacements for each refinement
for i, Nelem in enumerate(Plotting_element):
    
    ax1.plot(fem_displacement_domain_results[i], fem_displacement_results[i], label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical displacement solution
ax1.plot(analytical_domain, analytical_displacements, label='Analytical Displacement', color='black', linewidth=2)

ax1.legend()

# Plot Stresses
ax2.set_title('FEM Stresses vs Analytical Solution')

ax2.set_xlabel('Position along the bar (z) [mm]')

ax2.set_ylabel('Stress [MPa]')

ax2.grid(True)

# Plot FEM stresses for each refinement
for i, Nelem in enumerate(Plotting_element):
    
    x_values = fem_stress_domain_results[i]
    
    y_values = fem_stress_results[i]

    # Interpolate to match lengths if they are different
    if len(x_values) != len(y_values):
        # Create interpolation function
        interpolated_stress = np.interp(x_values, np.linspace(0, L, len(y_values)), y_values)

        # Plot interpolated stress
        ax2.plot(x_values, interpolated_stress, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')
    else:
        # Now both arrays are expected to have matching sizes
        ax2.plot(x_values, y_values, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical stress solution
ax2.plot(analytical_domain, analytical_stresses, label='Analytical Stress', color='black', linewidth=2)

ax2.legend()
plt.tight_layout()
plt.show()

#PLOT DISPLACMEENTS ERROR
plot_L2_error_vs_elements(Plotting_element, mse_disp_results)

plot_stress_error_vs_elements(Plotting_element, mse_results_stress)


# In[ ]:


import sympy as sp

import matplotlib.pyplot as plt

import numpy as np

#INPUT

A = 120 #mm^2
E = 70e3 #MPa
L = 500 #mm
a = 13 #N/mm
b = 0.13 #N/mm^2

#Calculating the length of each element

Length_element = L/Nelem;

#Boundary condition

BC = np.array([[1],[0]]);

#Symbolical z for evaluation of analytical relationship

z_symbol = sp.symbols('z');

#Introducing analytical solution ( I will get the values at the end of the bar)

analytical_end_displacement = 0.83829365

analytical_end_stress = (1 / A) * ((b*L**2)/2 + (a * L) - (a * L + (b*L**2)/2));

analytical_end_strain = analytical_end_stress / E;

# Function to calculate the relative error

def compute_relative_error_displacements(z_domain, analytical_disp, FEM_disp):
    
    #Create an empty list which will contain the values of the errors
    
    L_2_norm_error_disp = [];
    
    for node in range(int(len(z_domain))):
        
        error = np.linalg.norm(analytical_disp[node] - FEM_disp[node]);
        
        L_2_norm_error_disp.append(error);
        
        node += 1;
    
    return L_2_norm_error_disp

def compute_relative_error_stress(Nelem, analytical_stress, FEM_stress):
    
    #Create an empty list which will contain the values of the errors
    
    L_2_norm_error_stress = [];
    
    for element in range(Nelem):
        
        error = np.linalg.norm(analytical_stress[element] - FEM_stress[element]);
        
        L_2_norm_error_stress.append(error);
        
        element += 1;
    
    return L_2_norm_error_stress

#Define a library in order to plot the curves varying the number of element

Plotting_element = [1, 5, 10, 20, 30];

# Loop through different element sizes

for Nelem in Plotting_element:

    for element in range(Nelem):

        # Compute the FEM solution for the given number of elements

        Length_element = L / Nelem

        # Call your FEM functions here (assuming they return displacements):

        z_domain = mesh(Nelem, L)

        shape_functions = Assembly_shape_functions(Nelem, Length_element)

        derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)

        stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)

        connectivity_matrix = Assembly_connectivity_matrix(Nelem)

        global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)

        force_vector = Assembly_force_vector(Nelem, Length_element, a, b)

        new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

        FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)

        FEM_strain, FEM_stresses = Compute_stress_strain(FEM_displacements, Nelem, E)

    #Scale the displacements

    FEM_scaled_displacements, scaled_domain = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)
    
   #Compute the analytical displacements
    
    function_analytical_displacements_nodes = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a*z_symbol**2)/2 - (a*L + (b*L**2)/2)*z_symbol));

    numerical_function_analytical_disp = sp.lambdify(z_symbol, function_analytical_displacements_nodes);

    analytical_displacements = numerical_function_analytical_disp(z_domain);

    #Scale the stresses

    FEM_scaled_stresses, scaled_stresses_domain = scale_stresses(Nelem, shape_functions, FEM_displacements, Length_element, E)
    
    #Compute the analytical stresses

    function_analytical_stresses = -((1 / A) * ((b*z_symbol**2)/2 + (a * z_symbol) - (a * L + (b*L**2)/2)));

    numerical_analytical_stresses = sp.lambdify(z_symbol, function_analytical_stresses);

    analytical_stress = numerical_analytical_stresses(z_domain);

    #PLOT DISPLACEMENTS

    # Adjust the width and height of the plot

    plt.figure(figsize=(8, 4)) 

    #Plot FEM displacements along the length of the bar

    plt.plot(scaled_domain, FEM_scaled_displacements, marker='o', color='c', linestyle='solid')

    #Plot Analytical displacements along the length of the bar

    plt.plot(z_domain, analytical_displacements, color='r', linestyle=':')

    #Adjust the labels

    plt.xlabel("Length of the bar [mm]")

    plt.ylabel("Displacements")

    plt.title("Convergence Analysis (FEM vs Analytical Solution)")
    plt.grid(True)
    plt.show()

    #PLOT STRESS

    # Adjust the width and height of the plot

    plt.figure(figsize=(8, 4))

    #Plot FEM stresses along the length of the bar

    plt.plot(scaled_stresses_domain, FEM_scaled_stresses, marker='o', color='c', linestyle='solid')

    #Plot Analytical displacements along the length of the bar

    plt.plot(z_domain, analytical_stress, color='r', linestyle=':')

    #Adjust the labels

    plt.xlabel("Length of the bar [mm]")

    plt.ylabel("Stresses [MPa]")

    plt.title("Convergence Analysis (FEM vs Analytical Solution)")
    plt.grid(True)
    plt.show()




# In[19]:


len(scaled_stresses)


# In[ ]:


fem_displacements_results


# In[7]:


FEM_displacements


# In[11]:


int(len(analytical_domain))


# In[29]:


connectivity_matrix


# In[15]:


def stress_per_element(Nelem, shape_functions, displacements, Length_element, E):
    
    z_symbol = sp.symbols('z')  # Symbolic variable for position along the element
    
    element_stress = []
    element_numbers = []
    
    # Iterate through each element
    for elem in range(Nelem):
        # Get the nodal displacements for the current element (3 nodes per element)
        disp_node1 = displacements[2 * elem]  # Displacement at node 1
        disp_node2 = displacements[2 * elem + 1]  # Displacement at node 2
        disp_node3 = displacements[2 * elem + 2]  # Displacement at node 3
        
        # Combine displacements into a list
        displacements_elem = [disp_node1, disp_node2, disp_node3]
        
        # Compute the strain at the midpoint of the element (z = Length_element / 2)
        # Using the shape function derivatives at the midpoint to calculate strain
        midpoint_z = Length_element / 2
        strain_at_midpoint = 0
        
        # Compute the shape function derivatives for the midpoint
        for j in range(3):
            shape_func_derivative = sp.diff(shape_functions[f'Node {j + 1} Element {elem + 1}'], z_symbol)
            strain_at_midpoint += displacements_elem[j] * shape_func_derivative.subs(z_symbol, midpoint_z).evalf()
        
        # Compute stress using Hooke's law: stress = E * strain
        stress_at_midpoint = E * float(strain_at_midpoint)
        element_stress.append(stress_at_midpoint)
        
        # Store the element number (1-indexed)
        element_numbers.append(elem + 1)
    
    return element_stress, element_numbers

# Assuming you've defined 'shape_functions', 'displacements', and other necessary parameters
Nelem = 4  # Number of elements
Length_element = 10.0  # Length of each element
displacements = [0.0, 1.2, 2.5, 1.8, 3.4, 0.0]  # Example displacements
E = 200e9  # Young's modulus in Pascals

# Call the stress calculation function
element_stress, element_numbers = stress_per_element(Nelem, shape_functions, displacements, Length_element, E)

# Plot the results (stress against element number)
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))  # Set figure size
plt.plot(element_numbers, element_stress, marker='o', linestyle='-', color='b', label="Stress by Element")
plt.xlabel('Element Number')
plt.ylabel('Stress (Pa)')
plt.title('Stress Distribution by Element')
plt.legend()
plt.grid(True)
plt.xticks(range(1, Nelem + 1))  # Ensure x-axis ticks match the number of elements
plt.show()


# In[13]:


import matplotlib.pyplot as plt

# Define the range of element sizes for plotting
Plotting_element = [1, 5, 10, 20, 30]

# Initialize lists to store results for each refinement
fem_displacement_results = []
fem_domain_results = []

# Loop through different element sizes
for Nelem in Plotting_element:
    # Compute the FEM solution for the given number of elements
    Length_element = L / Nelem

    z_domain = mesh(Nelem, L)
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)
    force_vector = Assembly_force_vector(Nelem, Length_element, a, b)
    new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

    FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)
    FEM_scaled_displacements, scaled_domain = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)

    # Store FEM results for plotting later
    fem_displacement_results.append(FEM_scaled_displacements)
    fem_domain_results.append(scaled_domain)

# Compute the analytical displacements
function_analytical_displacements_nodes = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a*z_symbol**2)/2 - (a*L + (b*L**2)/2)*z_symbol))
numerical_function_analytical_disp = sp.lambdify(z_symbol, function_analytical_displacements_nodes)

# Generate analytical displacement values over the domain
analytical_domain = np.linspace(0, L, 500)  # More points for smoothness
analytical_displacements = numerical_function_analytical_disp(analytical_domain)

# Plot FEM results and analytical solution
plt.figure(figsize=(10, 6))

# Plot FEM displacements for each refinement
for i, Nelem in enumerate(Plotting_element):
    plt.plot(fem_domain_results[i], fem_displacement_results[i], label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical solution
plt.plot(analytical_domain, analytical_displacements, label='Analytical Solution', color='black', linewidth=2)

# Customize the plot
plt.title('FEM Displacements vs Analytical Solution')
plt.xlabel('Position along the bar (z)')
plt.ylabel('Displacement')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


# In[14]:


import numpy as np

import matplotlib.pyplot as plt

import sympy as sp

#INPUT

A = 120 #mm^2
E = 70e3 #MPa
L = 500 #mm
a = 13 #N/mm
b = 0.13 #N/mm^2

#Define a library in order to plot the curves varying the number of element

Plotting_element = [1, 2, 3, 4, 5];

#Boundary condition

BC = np.array([[1],[0]]);

#Symbolical z for evaluation of analytical relationship

z_symbol = sp.symbols('z');

# Initialize lists to store displacement and stress results for each refinement

fem_displacement_results = []

fem_displacement_domain_results = []

fem_stress_results = []

fem_stress_domain_results = []

# Loop through different element sizes

for Nelem in Plotting_element:
    
    # Compute the FEM solution for the given number of elements
    
    Length_element = L / Nelem

    z_domain = mesh(Nelem, L)
    
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    
    stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    
    connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    
    global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)
    
    force_vector = Assembly_force_vector(Nelem, Length_element, a, b)
    
    new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

    # Compute FEM displacements
    
    FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)

    # Scale displacements for plotting
    
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)
    
    # Compute FEM stresses
    
    FEM_strain, FEM_stresses = Compute_stress_strain(FEM_displacements, Nelem, E)
    
    #Scale stresses for plotting
    
    scaled_stresses, scaled_stress_domain = scale_stresses(Nelem, shape_functions, FEM_displacements, Length_element, E)
    
    # Store FEM displacement and stress results for plotting later
    
    fem_displacement_results.append(scaled_displacements)
    
    fem_displacement_domain_results.append(scaled_points)
    
    fem_stress_results.append(scaled_stresses)
    
    fem_stress_domain_results.append(scaled_stress_domain)  # z_domain corresponds to the element domain

# Analytical solution for displacements and stresses

# Displacement: u(z)

function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a*z_symbol**2)/2 - (a*L + (b*L**2)/2)*z_symbol))

numerical_function_analytical_displacement = sp.lambdify(z_symbol, function_analytical_displacement)

# Stress: σ(z) 

function_analytical_stress_nodes = -((1 / A) * ((b*z_symbol**2)/2 + (a * z_symbol) - (a * L + (b*L**2)/2)));

numerical_function_analytical_stress = sp.lambdify(z_symbol, function_analytical_stress_nodes)

# Generate analytical values over the domain

analytical_domain = np.linspace(0, L, 500)  # Fine grid for smooth plotting

analytical_displacements = numerical_function_analytical_displacement(analytical_domain)

analytical_stresses = numerical_function_analytical_stress(analytical_domain)

# Plot FEM and analytical solutions for displacements and stresses

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot Displacements
ax1.set_title('FEM Displacements vs Analytical Solution')

ax1.set_xlabel('Position along the bar (z) [mm]')

ax1.set_ylabel('Displacement [mm]')

ax1.grid(True)

# Plot FEM displacements for each refinement
for i, Nelem in enumerate(Plotting_element):
    
    ax1.plot(fem_displacement_domain_results[i], fem_displacement_results[i], label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical displacement solution
ax1.plot(analytical_domain, analytical_displacements, label='Analytical Displacement', color='black', linewidth=2)

ax1.legend()

# Plot Stresses
ax2.set_title('FEM Stresses vs Analytical Solution')

ax2.set_xlabel('Position along the bar (z) [mm]')

ax2.set_ylabel('Stress [MPa]')

ax2.grid(True)

# Plot FEM stresses for each refinement
for i, Nelem in enumerate(Plotting_element):
    
    x_values = fem_stress_domain_results[i]
    
    y_values = fem_stress_results[i]

    # Interpolate to match lengths if they are different
    if len(x_values) != len(y_values):
        # Create interpolation function
        interpolated_stress = np.interp(x_values, np.linspace(0, L, len(y_values)), y_values)

        # Plot interpolated stress
        ax2.plot(x_values, interpolated_stress, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')
    else:
        # Now both arrays are expected to have matching sizes
        ax2.plot(x_values, y_values, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical stress solution
ax2.plot(analytical_domain, analytical_stresses, label='Analytical Stress', color='black', linewidth=2)

ax2.legend()
plt.tight_layout()
plt.show()


# In[18]:


FEM_stresses


# In[11]:





# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Initialize lists to store displacement and stress results for each refinement
fem_displacement_results = []
fem_displacement_domain_results = []
fem_stress_results = []
fem_stress_domain_results = []

# Function to compute and scale stresses along the length of the bar
def scale_stresses(Nelem, shape_functions, displacements, Length_element, E):
    z_symbol = sp.symbols('z')  # Symbolic variable for position along the element
    
    element_stress = []
    element_centers = []  # to store midpoints
    
    # Iterate through each element
    for elem in range(Nelem):
        # Get the nodal displacements for the current element (3 nodes per element)
        disp_node1 = displacements[2 * elem]  # Displacement at node 1
        disp_node2 = displacements[2 * elem + 1]  # Displacement at node 2
        disp_node3 = displacements[2 * elem + 2]  # Displacement at node 3
        
        # Combine displacements into a list
        displacements_elem = [disp_node1, disp_node2, disp_node3]
        
        # Compute the strain at the midpoint of the element (z = Length_element / 2)
        midpoint_z = Length_element / 2
        
        strain_at_midpoint = 0
        
        # Compute the shape function derivatives for the midpoint
        for j in range(3):
            shape_func_derivative = sp.diff(shape_functions[f'Node {j + 1} Element {elem + 1}'], z_symbol)
            strain_at_midpoint += displacements_elem[j] * shape_func_derivative.subs(z_symbol, midpoint_z).evalf()
        
        # Compute stress using Hooke's law: stress = E * strain
        stress_at_midpoint = E * float(strain_at_midpoint)
        
        element_stress.append(stress_at_midpoint)
        element_centers.append(midpoint_z + elem * Length_element)  # Center position for the element
    
    return element_stress, element_centers

# Loop through different element sizes
for Nelem in Plotting_element:
    # Compute the FEM solution for the given number of elements
    Length_element = L / Nelem

    z_domain = mesh(Nelem, L)
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)
    force_vector = Assembly_force_vector(Nelem, Length_element, a, b)
    new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

    # Compute FEM displacements
    FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)

    # Scale displacements for plotting
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)

    # Compute FEM stresses
    FEM_stresses, stress_centers = scale_stresses(Nelem, shape_functions, FEM_displacements, Length_element, E)
    
    # Store FEM displacement and stress results for plotting later
    fem_displacement_results.append(scaled_displacements)
    fem_displacement_domain_results.append(scaled_points)
    fem_stress_results.append(FEM_stresses)
    fem_stress_domain_results.append(stress_centers)  # Use the stress centers for plotting

# Analytical solution for displacements and stresses
# Displacement: u(z) = a*z + (b*z^2)/2
function_analytical_displacement = a * z_symbol + (b * z_symbol**2) / 2
numerical_function_analytical_displacement = sp.lambdify(z_symbol, function_analytical_displacement)

# Stress: σ(z) = a + b*z
function_analytical_stress_nodes = a + b * z_symbol
numerical_function_analytical_stress = sp.lambdify(z_symbol, function_analytical_stress_nodes)

# Generate analytical values over the domain
analytical_domain = np.linspace(0, L, 500)  # Fine grid for smooth plotting
analytical_displacements = numerical_function_analytical_displacement(analytical_domain)
analytical_stresses = numerical_function_analytical_stress(analytical_domain)

# Plot FEM and analytical solutions for displacements and stresses
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot Displacements
ax1.set_title('FEM Displacements vs Analytical Solution')
ax1.set_xlabel('Position along the bar (z)')
ax1.set_ylabel('Displacement')
ax1.grid(True)

# Plot FEM displacements for each refinement
for i, Nelem in enumerate(Plotting_element):
    ax1.plot(fem_displacement_domain_results[i], fem_displacement_results[i], label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical displacement solution
ax1.plot(analytical_domain, analytical_displacements, label='Analytical Displacement', color='black', linewidth=2)

ax1.legend()

# Plot Stresses
ax2.set_title('FEM Stresses vs Analytical Solution')
ax2.set_xlabel('Position along the bar (z)')
ax2.set_ylabel('Stress')
ax2.grid(True)

# Plot FEM stresses for each refinement
for i, Nelem in enumerate(Plotting_element):
    x_values = fem_stress_domain_results[i]
    y_values = fem_stress_results[i]

    # Interpolate to match lengths if they are different
    if len(x_values) != len(y_values):
        # Create interpolation function
        interpolated_stress = np.interp(np.linspace(0, L, len(y_values)), x_values, y_values)

        # Plot interpolated stress
        ax2.plot(np.linspace(0, L, len(y_values)), interpolated_stress, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')
    else:
        # Now both arrays are expected to have matching sizes
        ax2.plot(x_values, y_values, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical stress solution
ax2.plot(analytical_domain, analytical_stresses, label='Analytical Stress', color='black', linewidth=2)

ax2.legend()
plt.tight_layout()
plt.show()


# In[9]:


# Compute the stresses using FEM
def compute_fem_stresses(Nelem, shape_functions, displacements, Length_element, E):
    z_symbol = sp.symbols('z')  # Symbolic variable for position along the element
    element_stress = []   # List to store stress values for each element
    element_centers = []  # List to store the center positions of each element

    # Iterate through each element
    for elem in range(Nelem):
        # Get the nodal displacements for the current element (3 nodes per element)
        disp_node1 = displacements[2 * elem]      # Displacement at node 1
        disp_node2 = displacements[2 * elem + 1]  # Displacement at node 2
        disp_node3 = displacements[2 * elem + 2]  # Displacement at node 3
        
        # Combine displacements into a list
        displacements_elem = [disp_node1, disp_node2, disp_node3]
        
        # Compute the strain at the midpoint of the element
        midpoint_z = (elem * Length_element) + (Length_element / 2)  # Center of the element
        
        strain_at_midpoint = 0
        
        # Compute the shape function derivatives for the midpoint
        for j in range(3):
            shape_func_derivative = sp.diff(shape_functions[f'Node {j + 1} Element {elem + 1}'], z_symbol)
            strain_at_midpoint += displacements_elem[j] * shape_func_derivative.subs(z_symbol, midpoint_z).evalf()
        
        # Compute stress using Hooke's law: stress = E * strain
        stress_at_midpoint = E * float(strain_at_midpoint)
        
        # Store the computed stress and center position
        element_stress.append(stress_at_midpoint)
        element_centers.append(midpoint_z)  # Store the center position for plotting
    
    return element_stress, element_centers

z_symbol = sp.symbols('z')

analytical_stress = -((1 / A) * ((b*z_symbol**2)/2 + (a * z_symbol) - (a * L + (b*L**2)/2)));



# Main function to run the simulation
def run_simulation(Nelem, L, a, b, E, A):
    Length_element = L / Nelem
    z_domain = mesh(Nelem, L)
    
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    
    # Assembly the stiffness matrix and force vector
    K_global = Assembly_global_stiffness_matrix(Nelem, Assembly_connectivity_matrix(Nelem), local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A))
    F_global = Assembly_force_vector(Nelem, Length_element, a, b)
    
    # Apply boundary conditions (example: fixed at the left end)
    BC = (np.array([1]), np.array([0]))  # Node 1 fixed at 0 displacement
    K_mod, F_mod = apply_boundary_conditions(K_global, BC, F_global)
    
    # Compute displacements
    displacements = Computation_of_displacements(K_mod, F_mod, Nelem, BC)
    
    # Compute FEM stresses
    fem_stresses, element_centers = compute_fem_stresses(Nelem, shape_functions, displacements, Length_element, E)
    
    # Compute analytical stresses for comparison
    analytical_values = [analytical_stress(a, b, L, z) for z in element_centers]
    
    # Plotting
    plt.figure(figsize=(10, 5))
    plt.plot(element_centers, fem_stresses, label='FEM Stress', marker='o')
    plt.plot(element_centers, analytical_values, label='Analytical Stress', linestyle='--')
    plt.xlabel('Position along the bar (z)')
    plt.ylabel('Stress')
    plt.title(f'Stress Distribution along the Bar (Nelem = {Nelem})')
    plt.legend()
    plt.grid()
    plt.show()

# Parameters
L = 10  # Length of the bar
a = 5   # Constant term for analytical solution
b = 1   # Linear term for analytical solution
E = 210e9  # Young's Modulus
A = 0.01  # Cross-sectional area

# Varying the number of elements
for Nelem in [2, 5, 10, 20]:
    run_simulation(Nelem, L, a, b, E, A)


# In[16]:


import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define the inputs
A = 120e-6  # mm^2 to m^2
E = 70e3  # MPa to Pa
L = 500e-3  # mm to m
a = 13  # N/mm to N/m
b = 0.13  # N/mm^2 to N/m^2

# Define the range of elements to plot
Plotting_element = [1, 2, 3, 4, 5]

# Boundary condition
BC = np.array([[1], [0]])

# Symbolical z for evaluation of analytical relationship
z_symbol = sp.symbols('z')

# Initialize lists to store displacement and stress results for each refinement
fem_displacement_results = []
fem_displacement_domain_results = []
fem_stress_results = []
fem_stress_domain_results = []

# Function to calculate L2 error
def calculate_L2_error(N_elements_range, analytical_func, fem_results, z_domain):
    L2_errors = []
    for i, Nelem in enumerate(N_elements_range):
        # Get the analytical values at the FEM domain points
        analytical_values = analytical_func(z_domain)

        # Compute L2 norm error
        error = np.sqrt(np.sum((analytical_values - fem_results[i]) ** 2) * (z_domain[1] - z_domain[0]))
        L2_errors.append(error)
    return L2_errors

# Loop through different element sizes
for Nelem in Plotting_element:
    # Compute the FEM solution for the given number of elements
    Length_element = L / Nelem
    z_domain = np.linspace(0, L, Nelem + 1)  # Mesh points along the bar

    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    stiffness_matrix = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    global_stiffness_matrix = Assembly_global_stiffness_matrix(Nelem, connectivity_matrix, stiffness_matrix)
    force_vector = Assembly_force_vector(Nelem, Length_element, a, b)
    
    new_global_stiffness, new_force_vector = apply_boundary_conditions(global_stiffness_matrix, BC, force_vector)

    # Compute FEM displacements
    FEM_displacements = Computation_of_displacements(new_global_stiffness, new_force_vector, Nelem, BC)

    # Scale displacements for plotting
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, FEM_displacements, Length_element)

    # Compute FEM stresses
    FEM_strain, FEM_stresses = Compute_stress_strain(FEM_displacements, Nelem, E)

    # Scale stresses for plotting
    scaled_stresses, scaled_stress_domain = scale_stresses(Nelem, shape_functions, FEM_displacements, Length_element, E)

    # Store FEM displacement and stress results for plotting later
    fem_displacement_results.append(scaled_displacements)
    fem_displacement_domain_results.append(scaled_points)
    fem_stress_results.append(scaled_stresses)
    fem_stress_domain_results.append(scaled_stress_domain)

# Analytical solution for displacements and stresses
function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a * z_symbol**2)/2 - (a * L + (b * L**2)/2) * z_symbol))
numerical_function_analytical_displacement = sp.lambdify(z_symbol, function_analytical_displacement)

function_analytical_stress_nodes = -((1 / A) * ((b * z_symbol**2) / 2 + (a * z_symbol) - (a * L + (b * L**2) / 2)))
numerical_function_analytical_stress = sp.lambdify(z_symbol, function_analytical_stress_nodes)

# Generate analytical values over the domain
analytical_domain = np.linspace(0, L, 500)  # Fine grid for smooth plotting
analytical_displacements = numerical_function_analytical_displacement(analytical_domain)
analytical_stresses = numerical_function_analytical_stress(analytical_domain)

# Plot FEM and analytical solutions for displacements and stresses
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot Displacements
ax1.set_title('FEM Displacements vs Analytical Solution')
ax1.set_xlabel('Position along the bar (z) [m]')
ax1.set_ylabel('Displacement [m]')
ax1.grid(True)

# Plot FEM displacements for each refinement
for i, Nelem in enumerate(Plotting_element):
    ax1.plot(fem_displacement_domain_results[i], fem_displacement_results[i], label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical displacement solution
ax1.plot(analytical_domain, analytical_displacements, label='Analytical Displacement', color='black', linewidth=2)
ax1.legend()

# Plot Stresses
ax2.set_title('FEM Stresses vs Analytical Solution')
ax2.set_xlabel('Position along the bar (z) [m]')
ax2.set_ylabel('Stress [Pa]')
ax2.grid(True)

# Plot FEM stresses for each refinement
for i, Nelem in enumerate(Plotting_element):
    x_values = fem_stress_domain_results[i]
    y_values = fem_stress_results[i]

    # Interpolate to match lengths if they are different
    if len(x_values) != len(y_values):
        # Create interpolation function
        interpolated_stress = np.interp(x_values, np.linspace(0, L, len(y_values)), y_values)
        # Plot interpolated stress
        ax2.plot(x_values, interpolated_stress, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')
    else:
        ax2.plot(x_values, y_values, label=f'FEM Nelem={Nelem}', marker='o', linestyle='--')

# Plot the analytical stress solution
ax2.plot(analytical_domain, analytical_stresses, label='Analytical Stress', color='black', linewidth=2)
ax2.legend()
plt.tight_layout()
plt.show()

# Calculate L2 errors
L2_errors = calculate_L2_error(Plotting_element, numerical_function_analytical_stress, fem_stress_results, analytical_domain)

# Plotting L2 errors
plt.figure(figsize=(10, 6))
plt.plot(Plotting_element, L2_errors, marker='o', color='blue', label='L2 Norm Error')
plt.xscale('linear')
plt.yscale('log')  # Log scale for better visibility of error reduction
plt.xlabel('Number of Elements')
plt.ylabel('L2 Norm Error (log scale)')
plt.title('L2 Norm Error Between Analytical and FEM Solution')
plt.legend()
plt.grid()
plt.show()


# In[48]:


def calculate_L2_norm_error_disp(Plotting_element, FEM_displacements, analytical_displacements_zdomain, L):
    
    # Create an array to store the L2 norm errors for each element count
    L2_error = np.zeros(len(Plotting_element))  # size based on number of element configurations

    # Iterate over the number of elements in Plotting_element
    for element_idx in range(len(Plotting_element)):
        
        # Skip if the displacements are not valid (i.e., 0 or non-array)
        if isinstance(analytical_displacements_zdomain[element_idx], (int, float)) or isinstance(FEM_displacements[element_idx], (int, float)):
            print(f"Skipping invalid displacement data for element index {element_idx}.")
            continue
        
        # Get the number of discretization points for the current element configuration
        N = len(analytical_displacements_zdomain[element_idx])

        # Discretization step size (depends on number of points in this configuration)
        delta_z = L / (N - 1)

        # Compute the L2 norm error using the volume integral approach for current element configuration
        L2_error[element_idx] = np.sqrt(
            np.sum((FEM_displacements[element_idx] - analytical_displacements_zdomain[element_idx])**2) * delta_z
        )
        
        # Example input checking code
        print(f"Type of FEM_displacements: {type(FEM_displacements)}")
        print(f"Type of analytical_displacements_zdomain: {type(analytical_displacements_zdomain)}")

    for i, arr in enumerate(analytical_displacements_zdomain):
        print(f"analytical_displacements_zdomain[{i}]: Type: {type(arr)}, Length: {len(arr)}, Values: {arr}")

    for i, arr in enumerate(FEM_displacements):
        print(f"FEM_displacements[{i}]: Type: {type(arr)}, Length: {len(arr)}, Values: {arr}")
    
    # Return the array of L2 norm errors for each element configuration
    return L2_error


# In[47]:


L2_error = calculate_L2_norm_error_disp(Plotting_element, FEM_displacements, analytical_displacements, L)


# In[10]:


L2_error = calculate_L2_error(z_domain, fem_displacement_results, analytical_displacements, Length_element)


# In[3]:


import numpy as np
import matplotlib.pyplot as plt

# Assuming fem_displacement_results and analytical_displacements are your results
fem_displacement_results = np.array(fem_displacement_results)  # FEM results
analytical_displacements = np.array(analytical_displacements)  # Analytical results

# Calculate L2 norm error for each node
l2_errors = np.sqrt(np.sum((fem_displacement_results - analytical_displacements) ** 2, axis=0))

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(l2_errors, marker='o', linestyle='-', color='b')
plt.title('L2 Norm Error between FEM and Analytical Displacements')
plt.xlabel('Node Index')
plt.ylabel('L2 Norm Error')
plt.grid()
plt.xticks(np.arange(len(l2_errors)))  # Set x-ticks to match node indices
plt.show()


# In[24]:


import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define parameters
L = 10  # Length of the bar
E = 210e9  # Young's Modulus in Pascals
A = 0.01  # Cross-sectional area in m^2
a = 1000  # Load parameter a
b = 2000  # Load parameter b

# Analytical displacement function (with z as the variable)
z_symbol = sp.symbols('z')
function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a * z_symbol**2)/2 - (a * L + (b * L**2) / 2) * z_symbol))

# Function to calculate L2-norm error
def calculate_l2_norm_error(fem_displacements, analytical_displacements):
    return np.sqrt(np.sum((fem_displacements - analytical_displacements) ** 2))

# Range of element counts
element_counts = np.arange(1, 21)  # From 1 to 20 elements
errors = []

for Nelem in element_counts:
    # Step 1: Mesh Generation
    z_domain = mesh(Nelem, L)
    
    # Step 2: Assemble Shape Functions and Derivatives
    Length_element = L / Nelem
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    
    # Step 3: Compute Local Stiffness Matrix
    local_stiffness_matrixes = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    
    # Step 4: Assemble Connectivity Matrix
    Connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    
    # Step 5: Assemble Global Stiffness Matrix
    Total_Element_stiffness_matrix = np.array(local_stiffness_matrixes)
    K_g = Assembly_global_stiffness_matrix(Nelem, Connectivity_matrix, Total_Element_stiffness_matrix)
    
    # Step 6: Assemble Force Vector
    F_global = Assembly_force_vector(Nelem, Length_element, a, b)
    
    # Step 7: Apply Boundary Conditions
    BC = np.array([[1], [0]])  # Example: First node fixed (u = 0)
    K_mod, F_mod = apply_boundary_conditions(K_g, BC, F_global)
    
    # Step 8: Compute Displacements
    displacements = Computation_of_displacements(K_mod, F_mod, Nelem, BC)
    
    # Step 9: Scale Displacements
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, displacements, Length_element)
    
    # Step 10: Analytical Displacements at FEM nodes
    analytical_displacements = np.array([function_analytical_displacement.subs(z_symbol, z).evalf() for z in z_domain])
    
    # Step 11: Calculate L2-norm error
    error = calculate_l2_norm_error(displacements, analytical_displacements)
    errors.append(error)

# Plotting the L2-norm error
plt.figure(figsize=(10, 6))
plt.plot(element_counts, errors, marker='o')
plt.title('L2-norm Error vs Number of Elements')
plt.xlabel('Number of Elements')
plt.ylabel('L2-norm Error')
plt.grid()
plt.yscale('log')  # Optional: use logarithmic scale for better visualization
plt.xticks(element_counts)
plt.show()


# In[23]:


analytical_displacements


# In[29]:


import numpy as np
import matplotlib.pyplot as plt

# Function to calculate L2 norm error
def compute_L2_norm_error(fem_values, analytical_values, Length_element):
    return np.sqrt(np.sum((fem_values - analytical_values) ** 2) * Length_element)

def calculate_L2_error(Nelem, fem_values, analytical_values, Length_element):
    """
    Computes the L2 norm error for given FEM and analytical values.

    Parameters:
    - Nelem: Number of elements used in the FEM analysis.
    - fem_values: The displacements computed from the FEM.
    - analytical_values: The corresponding analytical displacement values.
    - Length_element: The length of each element.

    Returns:
    - L2_error: The computed L2 norm error.
    """
    # Ensure both arrays are numpy arrays
    fem_values = np.array(fem_values)
    analytical_values = np.array(analytical_values)

    # Calculate the L2 error
    L2_error = compute_L2_norm_error(fem_values, analytical_values, Length_element)
    
    return L2_error

def plot_L2_error_vs_elements(Nelem_list, fem_displacements_list, analytical_displacements, Length):
    """
    Plots the L2 error vs number of elements for FEM and analytical solutions.

    Parameters:
    - Nelem_list: List of number of elements.
    - fem_displacements_list: List of FEM displacement arrays (one for each Nelem).
    - analytical_displacements: The analytical displacements array.
    - Length: Length of the bar (domain length).
    """
    # List to store L2 error values
    L2_errors = []

    # Loop over each element count (Nelem)
    for Nelem, fem_values in zip(Nelem_list, fem_displacements_list):
        
        # Interpolate analytical displacements to match fem_values size
        
        analytical_values = np.interp(np.linspace(0, Length, len(fem_values)),
                                      np.linspace(0, Length, len(analytical_displacements)),
                                      analytical_displacements)

        # Compute the length of each element
        Length_element = Length / Nelem
        
        # Compute the L2 norm error
        L2_error = calculate_L2_error(Nelem, fem_values, analytical_values, Length_element)
        
        # Store the L2 error
        L2_errors.append(L2_error)

    # Plot the L2 error vs number of elements
    plt.figure(figsize=(8, 6))
    plt.loglog(Nelem_list, L2_errors, marker='o', linestyle='-', color='b', label='L2 Error')
    plt.xlabel('Number of Elements')
    plt.ylabel('L2 Error')
    plt.title('L2 Error vs Number of Elements (FEM vs Analytical)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()


# In[30]:


plot_L2_error_vs_elements(Plotting_element, FEM_displacements, analytical_displacements_zdomain, L)


# In[28]:


len(FEM_displacements)


# In[31]:


import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define parameters
L = 10  # Length of the bar
E = 210e9  # Young's Modulus in Pascals
A = 0.01  # Cross-sectional area in m^2
a = 1000  # Load parameter a
b = 2000  # Load parameter b

# Analytical displacement function (with z as the variable)
z_symbol = sp.symbols('z')
function_analytical_displacement = -((1 / (E * A)) * ((b * z_symbol**3)/6 + (a * z_symbol**2)/2 - (a * L + (b * L**2) / 2) * z_symbol))

# Function to calculate L2-norm error
def calculate_l2_norm_error(fem_displacements, analytical_displacements):
    return np.sqrt(np.sum((fem_displacements - analytical_displacements) ** 2))

# Range of element counts
element_counts = np.arange(1, 21)  # From 1 to 20 elements
errors = []

for Nelem in element_counts:
    # Step 1: Mesh Generation
    z_domain = mesh(Nelem, L)
    
    # Step 2: Assemble Shape Functions and Derivatives
    Length_element = L / Nelem
    shape_functions = Assembly_shape_functions(Nelem, Length_element)
    derivatives_shape_functions = compute_derivatives_of_shape_functions(shape_functions)
    
    # Step 3: Compute Local Stiffness Matrix
    local_stiffness_matrixes = local_stiffness_matrix(derivatives_shape_functions, Nelem, Length_element, E, A)
    
    # Step 4: Assemble Connectivity Matrix
    Connectivity_matrix = Assembly_connectivity_matrix(Nelem)
    
    # Step 5: Assemble Global Stiffness Matrix
    Total_Element_stiffness_matrix = np.array(local_stiffness_matrixes)
    K_g = Assembly_global_stiffness_matrix(Nelem, Connectivity_matrix, Total_Element_stiffness_matrix)
    
    # Step 6: Assemble Force Vector
    F_global = Assembly_force_vector(Nelem, Length_element, a, b)
    
    # Step 7: Apply Boundary Conditions
    BC = np.array([[1], [0]])  # Example: First node fixed (u = 0)
    K_mod, F_mod = apply_boundary_conditions(K_g, BC, F_global)
    
    # Step 8: Compute Displacements
    displacements = Computation_of_displacements(K_mod, F_mod, Nelem, BC)
    
    # Step 9: Scale Displacements
    scaled_displacements, scaled_points = scaling_displacements(z_domain, Nelem, shape_functions, displacements, Length_element)
    
    # Step 10: Analytical Displacements at FEM nodes
    analytical_displacements = np.array([function_analytical_displacement.subs(z_symbol, z).evalf() for z in z_domain])
    
    # Step 11: Calculate L2-norm error
    error = calculate_l2_norm_error(scaled_displacements, analytical_displacements)
    errors.append(error)

# Plotting the L2-norm error
plt.figure(figsize=(10, 6))
plt.plot(element_counts, errors, marker='o')
plt.title('L2-norm Error vs Number of Elements')
plt.xlabel('Number of Elements')
plt.ylabel('L2-norm Error')
plt.grid()
plt.yscale('log')  # Optional: use logarithmic scale for better visualization
plt.xticks(element_counts)
plt.show()


# In[ ]:




