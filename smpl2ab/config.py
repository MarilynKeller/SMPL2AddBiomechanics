import os

package_directory = os.path.dirname(os.path.abspath(__file__))

smpl_folder = os.path.join(package_directory, '../models')
osim_model_path = os.path.join(package_directory, '../models/bsm/bsm.osim')
osim_sample_motion = os.path.join(package_directory, '../models/bsm/sample_motion/01/01_01_ik.mot')
osim_geometry_folder = os.path.join(package_directory, '../models/bsm/Geometry')

# Marker sets
bsm_markers_on_smpl_path = os.path.join(package_directory, "data/bsm_markers.yaml")
bsm_markers_on_smplx_path = os.path.join(package_directory, "data/bsm_markers_smplx.yaml")

mesh_vertices_path = os.path.join(package_directory, 'measurements/smpl_measurement_vertices.yaml' )    
mesh_vertices_smplx_path = os.path.join(package_directory, 'measurements/smpl_measurement_vertices_smplx.yaml' )    
bsm_osso_segmentation = os.path.join(package_directory, "data/OSSO_osim_bone_groups_to_vertices_v4.pkl") #Dict that stores the segmentation of OSSO into rajagopal bones. This was created with Blender.
osso_inf_template = os.path.join(package_directory, "data/osso_inference_topo.ply") #Template mesh generated by the OSSO inference
smpl_markers_folder = '/tmp/'
osso_rj_unposed = os.path.join(package_directory, "data/osso_rajagopal_unposed_v2.obj")
osim_unposed = os.path.join(package_directory, "data/rajagopal_unposed.obj")