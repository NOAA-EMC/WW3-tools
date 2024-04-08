import gzip
import tarfile
import os

def unzip_and_untar(gz_file, output_dir):
    # Extracting the filename without extension
    file_name = os.path.splitext(os.path.basename(gz_file))[0]
    print(file_name)
    # Get the base name of the extracted folder
    base_name = file_name.rsplit('.', 1)[0]  # Get everything before the last dot

    # Setting the output path for the extracted file
    output_file = os.path.join(output_dir, file_name)

    # Unzipping the .gz file
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            f_out.write(f_in.read())

    # Untaring the extracted file
    with tarfile.open(output_file, 'r') as tar:
        tar.extractall(output_dir)

    # Get the base name of the extracted folder
    extracted_folder = os.path.join(output_dir, base_name)

    print(extracted_folder)
    # Creating a list of the extracted files
    list_file = os.path.join(output_dir, f'{base_name}_contents.txt')
    with open(list_file, 'w') as f:
        for root, dirs, files in os.walk(extracted_folder):
            for file in files:
                f.write(os.path.join(root, file) + '\n')

    # Creating a list of ids between two dots in the filenames
    id_file = os.path.join(output_dir, f'{base_name}_id.txt')
    with open(id_file, 'w') as f:
        for root, dirs, files in os.walk(extracted_folder):
            for file in files:
                file_parts = file.split('.')
                if len(file_parts) >= 3:
                    id_between_dots = file_parts[1]
                    f.write(id_between_dots + '\n')

# Define the input .gz file and the output directory
gz_file = 'multi_1.t00z.spec_tar.gz'
output_dir = '.'  # You can specify any output directory here

# Call the function to unzip and untar the file
unzip_and_untar(gz_file, output_dir)

