"""
To shift the mesh from -180:180 (lon) to 0:360

input_file_path: -180:180 (lon) mesh in GMSH fomrat
output_fil_pat: 0:360 (lon) mesh in GMSH format 
"""

# Author: Ali Salimi-Tarazouj


input_file_path = './9km_nobc.ww3'
output_file_path = './modified_9km_nobc.ww3'

# Initialize lists to store node data
node_data = []

# Read the input mesh file, extract and adjust node data
with open(input_file_path, 'r') as input_file:
    lines = input_file.readlines()
    node_section = False

    try:
        start_node_index = lines.index('$Nodes\n') + 1
        end_node_index = lines.index('$EndNodes\n')

        for line in lines[start_node_index:end_node_index]:
            line = line.strip()  # Remove leading/trailing whitespace
            parts = line.split()

            # Check if there are enough values in the line (at least 4)
            if len(parts) >= 4:
                x, y, z = map(float, parts[1:4])
                if x < 0:
                    x += 360
                node_data.append((x, y, z))

        # Write the modified mesh to the output file
        with open(output_file_path, 'w') as output_file:
            output_file.write(''.join(lines[:start_node_index]))
            output_file.write(f'{len(node_data)}\n')
            for i, (x, y, z) in enumerate(node_data, start=1):
                output_file.write(f'{i} {x:.5f} {y:.5f} {z:.5f}\n')
            output_file.write(''.join(lines[end_node_index:]))
    
    except ValueError:
        print("Error: '$Nodes' or '$EndNodes' not found in the input file.")

