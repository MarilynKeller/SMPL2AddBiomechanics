import os


package_directory = os.path.dirname(os.path.abspath(__file__))

smpl_folder = os.path.join(package_directory, '../models')
bsm_model_path = os.path.join(package_directory, '../models/bsm.osim')

# Do not edit the following lines

mesh_vertices_path = os.path.join(package_directory, 'measurements/smpl_measurement_vertices.yaml' )    
bsm_osso_segmentation = "/is/cluster/work/mkeller2/Data/TML/OSSO/OSSO_osim_bone_groups_to_vertices_v4.pkl" #Dict that stores the segmentation of OSSO into rajagopal bones. This was created with Blender.