import numpy as np

import math as mt

def calcUSSR(E,BC,F,Con,NodePos,Area):
    
    import numpy as np

    import math as mt

    def Computation_Lengths(Con, NodePos):

        Nelem = np.shape(Con)[0];

        Difference_y = np.zeros(Nelem);

        Difference_x = np.zeros(Nelem);

        Length = np.zeros(Nelem);

        #Take the x-coordinates and y-coordinates for each node

        y_coordinates = NodePos[:, 1];

        x_coordinates = NodePos[:, 0];

        start_point_indices = Con[:, 0];

        end_point_indices = Con[:, 1];

        #Computation of the differences between x-coordinates and y-coordinates for each element

        for j in range(Nelem):

            start_index = start_point_indices[j];

            end_index = end_point_indices[j];

            Difference_y[j] = y_coordinates[int(end_index - 1)] - y_coordinates[int(start_index - 1)];

            Difference_x[j] = x_coordinates[int(end_index - 1)] - x_coordinates[int(start_index - 1)];

            j += 1;

        #Computation of the length of each element of the system

        for i in range(len(Difference_y)):

            Length[i] = np.sqrt((Difference_x[i]**2) + (Difference_y[i]**2));

            i += 1;

        return Length


    def Computation_of_local_stiffness(Con, Area, Length, E):

        Nelem = np.shape(Con)[0];

        Global_Stiffness_local = np.zeros([int(Nelem), int(Nelem)]);

        Local_element_stiffness = np.zeros(Nelem);

        #Computation of the local stiffness matrix of each element

        for i in range(Nelem): 

            Local_element_stiffness[i] = E*Area[i]/Length[i];

        return Local_element_stiffness


    def Assemble_the_element_stiffness_matrixes(Local_element_stiffness, Con):

        Nelem = np.shape(Con)[0];

        #Put all the element stiffness matrixes in one dictionary so it will be easier to recall them later on

        Total_Element_stiffness_matrix = {};

        #Keys ---- > 'Element 1,2,..'

        keys = [];

        #Values -----> Element stiffness matrixes

        values = [];

        temp = 0;

        for rows in range(int(np.shape(Con)[0])):

            #Initialization of the element stiffness matrix

            Element_stiffness_matrix = np.zeros([4, 4]);

            for j in range(4):

                for z in range(4):

                    #1 CASE: The element in the matrix is found in the x-positions of the matrix where i = j:

                    if (j == 0 or j == range(4)[-2]) and (z == 0 or z == range(4)[-2]) and ( j == z ):

                        Element_stiffness_matrix[j][z] += Local_element_stiffness[temp];

                        z += 1

                    #2 CASE: The element in the matrix is found in the x-positions of the matrix where i != j

                    elif (j == 0 or j == range(4)[-2]) and (z == 0 or z == range(4)[-2]) and ( j != z ):

                        Element_stiffness_matrix[j][z] -= Local_element_stiffness[temp];

                        z += 1;

                    #3 CASE: The element in the matrix is found in the y-columns (second and fourth)

                    elif (j == 1 or j == range(4)[-1]) and (z == 1 or z == range(4)[-1]):

                        Element_stiffness_matrix[j][z] = 0;

                        z += 1;

                j += 1;

            temp += 1;

            #Add the stiffness matrix to Values

            values.append(Element_stiffness_matrix);

            #Add 'Element (rows + 1)' to Keys

            keys.append(f'Element {rows + 1}');

        #Assemble the dictionary

        for key, value in zip(keys,values):

            Total_Element_stiffness_matrix[key] = value;      

        return Total_Element_stiffness_matrix


    def Rotation_of_element_stiffness_matrixes(Total_Element_stiffness_matrix, Con, NodePos):              

        Nelem = np.shape(Con)[0];

        #Get the indices of the starting point of all elements

        start_point_indices = Con[:,0];

        #Get the endices of the end points of all elements

        end_point_indices = Con[:,1];

        #Initialization of the rotational matrix of all elements

        Total_rotational_matrix = {};

        #Initialization of the rotated stiffness matrixes for all elements

        Total_rotated_matrixes = {};

        keys = [];

        rotations = [];

        rotated_matrixes = [];

        theta = np.zeros(Nelem)

        #Inizialization of the angles index( for the calculation of the rotational matrix )

        temp = 0;

        #Computation of the angles of each element

        for i in range(int(np.shape(Con)[0])):

            for j in range(int(np.shape(Con)[0])):

                First_node = NodePos[int(start_point_indices[i] - 1)];

                Second_node = NodePos[int(end_point_indices[i] - 1)];

                dy = (Second_node[1] - First_node[1]);

                dx = (Second_node[0] - First_node[0]);

                if dx == 0:

                    theta[i] = 0;

                else:

                    theta[i] = mt.degrees(mt.atan(dy/dx));

                if theta[i] < 0:

                    theta[i] = 180 + theta[i];

                j += 1;

            keys.append(f'Element {i + 1}');

            #Initialization of an elementary rotational matrix

            Element_Rotational_matrix = [[mt.cos(mt.radians(theta[i])), -(mt.sin(mt.radians(theta[i]))), 0, 0], [mt.sin(mt.radians(theta[i])), mt.cos(mt.radians(theta[i])), 0, 0], [0, 0, mt.cos(mt.radians(theta[i])), -(mt.sin(mt.radians(theta[i])))], [0, 0, mt.sin(mt.radians(theta[i])), mt.cos(mt.radians(theta[i]))]];

            #Inverting the rotational matrix

            inverted_rotational_matrix = np.linalg.inv(Element_Rotational_matrix);

            #    K_rotated = T * K * [T]-1                  Recalling the element stiffness matrix

            rotated_matrix = Element_Rotational_matrix@Total_Element_stiffness_matrix[f'Element {i + 1}']@inverted_rotational_matrix;


            rotations.append(Element_Rotational_matrix);

            rotated_matrixes.append(rotated_matrix);

            i += 1;

        for key, value in zip(keys, rotations):

            Total_rotational_matrix[key] = value;

        for key, matrix in zip(keys, rotated_matrixes):

            Total_rotated_matrixes[key] = matrix;

        return Total_rotational_matrix, Total_rotated_matrixes

    def Assembly_global_stiffness_matrix(Total_Element_stiffness_matrix, Con, NodePos):

            Nelem = np.shape(Con)[0];

            #Initialization of the global stiffness matrix

            K_g = np.zeros([2*(NodePos.shape[0]), 2*(NodePos.shape[0])]);

            for index in range(Nelem):

                #Recalling the stiffness matrix

                stiffness_matrix = Total_Element_stiffness_matrix[f'Element {index + 1}'];

                #Get the start point of the element

                start_point = Con[index, 0];

                #Get the end point of the element

                end_point = Con[index, 1];

                #Adjustement of the indexes (Relationships between the nodes of the element stiffness matrix and those of the global one)

                new_indexes = [(2*start_point) - 2, (2*start_point) - 1, (2*end_point) - 2, (2*end_point) - 1];

                for i in range(4):

                    for j in range(4):

                        K_g[new_indexes[i]] [new_indexes[j]] += stiffness_matrix[i][j];

                        j += 1;

                    i += 1;

                index += 1

            return K_g

    def apply_boundary_conditions(K_g, BC, F, NodePos):
        K_mod = K_g.copy()  # Copy the global stiffness matrix
        F_mod = F.copy()

        nodes_with_disp = BC[0] - 1  # Adjust to zero-based indexing
        displacements = BC[1]

        sorted_indices = np.argsort(nodes_with_disp)[::-1]  # Sort nodes for deletion
        nodes_with_disp = nodes_with_disp[sorted_indices]
        displacements = displacements[sorted_indices]

        colrem = np.zeros(len(F_mod))  # Column reminder for stiffness

        for i, node in enumerate(nodes_with_disp):
            disp = displacements[i]
            if disp == 0:  # Nodes with zero displacement
                K_mod = np.delete(K_mod, node, axis=0)
                K_mod = np.delete(K_mod, node, axis=1)
                F_mod = np.delete(F_mod, node)
                colrem = np.delete(colrem, node)
            else:  # Nodes with non-zero displacement
                if node < K_mod.shape[1]:
                    colrem[:len(K_mod)] -= K_mod[:, node] * disp
                K_mod = np.delete(K_mod, node, axis=0)
                K_mod = np.delete(K_mod, node, axis=1)
                F_mod = np.delete(F_mod, node)
                colrem = np.delete(colrem, node)

        result = F_mod + colrem[:len(F_mod)]

        if K_mod.size > 0:
            Finaldisplacements = np.linalg.solve(K_mod, result)

        else:
            Finaldisplacements = np.array([])  # No displacements if all nodes constrained

        total_displacements = np.zeros(2 * NodePos.shape[0])  # Total displacements vector

        for i, node in enumerate(BC[0]):
            total_displacements[node - 1] = BC[1][i]

        free_nodes = [i for i in range(2 * NodePos.shape[0]) if i not in (BC[0] - 1)]
        for i, node in enumerate(free_nodes):
            total_displacements[node] = Finaldisplacements[i]

        return K_mod, result, total_displacements

    def calcStrain(Con, total_displacements, NodePos, Total_Rotational_matrix, Length):

        Nelem = np.shape(Con)[0];

        Total_strains = np.zeros(Nelem)

        strains_x = np.zeros(Nelem);

        strains_y = np.zeros(Nelem);

        for index in range(Nelem):

            start_point = Con[index][0];

            end_point = Con[index][1];

            rotated_considered_disp = np.zeros(4);

            Total_considered_disp = np.zeros(4);

            considered_disp = np.zeros(4);

            new_index = [(2*start_point)-2, (2*start_point)-1, (2*end_point)-2, (2*end_point)-1];

            #Get the displacements for each node of the element in the connectivity

            for disp in range(4):

                rotated_considered_disp[disp] += total_displacements[int(new_index[disp])];

                disp += 1;

            #Rotation of the displacements

            Rotation_matrix = Total_Rotational_matrix[f'Element {index + 1}'];

            Inverted_rotation_matrix = np.linalg.inv(Rotation_matrix);

            Total_considered_disp = Inverted_rotation_matrix@rotated_considered_disp;

            #Consider only the x-components

            considered_disp[0] = Total_considered_disp[0];

            considered_disp[1] = Total_considered_disp[2];

            #Computation of the strains

            Total_strains[index] = (considered_disp[1] - considered_disp[0])/Length[index]


        return Total_strains 

    def calcStress(Nelem, Total_strains, E):

        stress = np.zeros(Nelem)

        for i in range(len(stress)):

            stress[i] = E*Total_strains[i];

            i += 1;

        return stress

    def calc_Reaction_Forces(BC, K_g, F, Con, total_displacements, NodePos):

        Nelem = np.shape(Con)[0];

        reaction_nodes = BC[:][0];

        reaction_forces = np.zeros(int((2*NodePos.shape[0])));

        reaction_forces = ((K_g@total_displacements) - F);

        #Eliminating the small errors

        for i in range(len(total_displacements)):

            if total_displacements[i] != 0:

                reaction_forces[i] = 0;

            i += 1;


        return reaction_forces

    Nelem = np.shape(Con)[0];

    Length = Computation_Lengths(Con, NodePos)

    Local_element_stiffness = Computation_of_local_stiffness(Con, Area, Length, E)

    Local_element_stiffness_matrix = Assemble_the_element_stiffness_matrixes(Local_element_stiffness, Con)

    Rotational_matrixes, Rotated_element_stiffness_matrixes = Rotation_of_element_stiffness_matrixes(Local_element_stiffness_matrix, Con, NodePos)

    K_g = Assembly_global_stiffness_matrix(Rotated_element_stiffness_matrixes, Con, NodePos)

    K_g_applied, Force_vector_new, displacements = apply_boundary_conditions(K_g, BC, F, NodePos)

    Reaction_forces = calc_Reaction_Forces(BC, K_g, F, Con, displacements, NodePos);

    Strains = calcStrain(Con, displacements, Nelem, Rotational_matrixes, Length);

    Stresses = calcStress(Nelem, Strains, E);

    return displacements, Strains, Stresses, Reaction_forces

