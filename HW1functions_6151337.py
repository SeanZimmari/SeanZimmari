import numpy as np


def calcDisp(L, Nelem, F_ext, E, BC):
   
    def CalcLelem(L, Nelem):  # Calculate element length
        return L / Nelem

    LE = CalcLelem(L, Nelem)

   
    def Area(L, Nelem): #  Calculate areas for each element
        nodepos = np.linspace(0, L, Nelem + 1)
        w1 = 50
        w2 = 25
        t = 3.125
        areas = (w1 + ((w2 - w1) / L) * nodepos) * t
        return areas[1:]  # Element areas for Nelem without first node 

    elemArea = Area(L, Nelem)

  
    def keq(areas, E, Lelem):   # Calculate the equivalent stiffness for each element
        scalar = E / Lelem
        keq_values = areas * scalar  # Element stiffness values
        return keq_values

    keq_values = keq(elemArea, E, LE)

    
    def glob_stiff_matrix(Nelem):   # Construction of the global stiffness matrix with connectivity matrix
        
        local_k_matrix = [np.array([[k, -k], [-k, k]]) for k in keq_values] #Local stiff matrices
        CM = np.array([[i, i+1] for i in range(1, Nelem + 1)]) # Connectivity matrix
        K_G = np.zeros((Nelem + 1, Nelem + 1))

        for i in range(Nelem): # Adjust indexing
            raw = CM[i]
            newraw = [n - 1 for n in raw]  

            for a in range(2): #Connecting each node to each local stiff matrix and filling the global stiffness matrix
                for b in range(2):
                    pos1 = newraw[a]
                    pos2 = newraw[b]
                    K_G[pos1, pos2] += local_k_matrix[i][a][b] # Select local stifness matrix and elements in it 
        return K_G

    K_G = glob_stiff_matrix(Nelem)

    
    def apply_boundary_conditions(K_G, BC, F_ext): # Setting buondary condition for zero and non-zero displacements
    
        F_global = np.zeros(Nelem + 1)
        F_global[-1] = F_ext  # Apply external force at the last node of force vector

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

  
    K_mod, result = apply_boundary_conditions(K_G, BC, F_ext) # BC application to solve for the unknown displacements

    if K_mod.size > 0:
        Finaldisplacements = np.linalg.solve(K_mod, result)
    else:
        Finaldisplacements = np.array([])  # If BC applied to all the nodes then no displacements to find

   
    total_displacements = np.zeros(Nelem + 1) # Recreating the original dimension vector (Nelem+1) to store all displacements including free nodes and boundary conditions

   
    for i, node in enumerate(BC[0]):# Add to the total displacements vector the displacements of the nodes with BC by re-indexing
        total_displacements[node - 1] = BC[1][i]

    
    free_nodes = [i for i in range(Nelem + 1) if i not in (BC[0] - 1)] # Free nodes indexing from BC, then add to the total displacements vector the found free displacements 
    for i, node in enumerate(free_nodes):
        total_displacements[node] = Finaldisplacements[i]

    strains = np.zeros(Nelem) 
    stresses = np.zeros(Nelem)

    for i in range(Nelem):# Calculate stress and strain for each node and element
        strains[i] = (total_displacements[i + 1] - total_displacements[i]) / LE
        stresses[i] = E * strains[i]

    return total_displacements, strains, stresses