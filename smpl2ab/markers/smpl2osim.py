import pickle
import numpy as np

import trimesh
from markers.marker_transfer import SSMarkerTransfer
import nimblephysics as nimble
import config as cg
from smpl2ab.markers.osim_editor import OsimEditor
from aitviewer.renderables.osim import OSIMSequence
from aitviewer.viewer import Viewer

class Smpl2osim:
    
    def __init__(self, smpl_marker_dict:dict, 
                 osim:nimble.biomechanics.OpenSimFile, 
                 osso_segmentation:dict, 
                 osim_model_path:str, 
                 marker_set_name:str, 
                 rigging_method:str, 
                 mapping_method='per_part') -> None:
        """_summary_

        Args:
            smpl_marker_dict (dict): Dictionary of markers in the SMPL model. {'markerlabel1': vertex_index1, 'markerlabel2': vertex_index2, ...}
            osim_model (_type_): OpenSim skeleton model with .osim extension
            osso_segmentation (dict): Segmentation of the mean OSSO body, matching the osim model bones.
            osim_model_path (str): Path to the osim model
            rigging_method (str): Method to to determine, for each marker to which osim bone to rig it. Possible: closest_rj_bone, heavy_smpl_bone, gt
        """
        
        self.smpl_marker_dict = smpl_marker_dict
        self.osso_segmentation = osso_segmentation
        self.osim = osim
        self.marker_set_name = marker_set_name
        self.osim_model_path = osim_model_path
        
        self.rigging_method = rigging_method
        self.mapping_method = mapping_method

        osim_node_names = [n.getName() for n in osim.skeleton.getBodyNodes()]
        osim_node_names.sort()
        
        # Check that the segmentation is valid
        for k in osso_segmentation.keys():
            if k not in osim_node_names:
                raise ValueError(f'Segmentation of OSSO contains {k}, which is not a bone in the osim model. Osim bones are: {[bone for bone in osim_node_names]}')                
                
        print(f'Marker numbers = {len(self.smpl_marker_dict)}')
        print(f'Osim model bones = {osim_node_names}')
        
        
    @classmethod
    def from_files(cls, smpl_marker_dict:dict, osim_model_path:str, osso_segmentation:str, **kwargs) -> 'Smpl2osim':
    
        osim_model = nimble.biomechanics.OpenSimParser.parseOsim(osim_model_path)
        
        return cls(smpl_marker_dict, osim_model, osso_segmentation, osim_model_path, **kwargs)
   
    def generate_osim(self, smpl_mesh_path, osso_mesh_path, output_osim_path, display=False):

        # Todo, pick a per gender mean body or target body
        # Load a template skin
        smpl_trimesh = trimesh.load_mesh(smpl_mesh_path, process=False)
        # Apply the markers on the skin
        skin_indices = np.hstack(self.smpl_marker_dict.values())
        markers_loc = smpl_trimesh.vertices[skin_indices]
        markers_label = [k for k in self.smpl_marker_dict.keys()]   
        
        # Load the corresponding OSSO skeleton
        osso_p = trimesh.load_mesh(osso_mesh_path, process=False)
        osso_rj = trimesh.load_mesh(cg.osso_rj_unposed, process=False, maintain_order=True, skip_materials=True)

        # Compute the relative position of markers and rigging
        mt = SSMarkerTransfer(markers_loc, osso_p, osso_rj, markers_label, skin_indices, rigging_method=self.rigging_method, mapping_method=self.mapping_method)
        mt.export_correspondances(self.marker_set_name)
        if display:
            print('Displaying the marker correspondances')
            mt.visualize(smpl_trimesh)
        
        # Edit the osim with computed markers and save the resulting osim model
        oe = OsimEditor(self.osim_model_path, self.marker_set_name)
        oe.export_osim(output_osim_path)
        if display:
            v = visualize_osim(output_osim_path, self.marker_set_name)
            v.run()
            
            

def visualize_osim(osim_path, marker_set):


    v = Viewer()
    osim_model_results = OSIMSequence.a_pose(osim_path=osim_path, name='result_osim', color_markers_per_part=True, color_markers_per_index=False, color_skeleton_per_part=True)
    v.scene.add(osim_model_results)


    # 
    # # import ipdb; ipdb.set_trace()
    # c = (255/255, 85/255, 255/255, 1)
    # from aitviewer.utils import vertex_colors_from_weights
    # # markers_pos_abs = markers_pos_abs[:100,:]
    # colors = vertex_colors_from_weights(weights=range(len(markers_pos_abs)), scale_to_range_1=True, alpha=1)[np.newaxis, :, :]
    # markers_abs_pc = PointClouds(markers_pos_abs, colors=colors, position = [0,0,1], name='markers_abs', point_size=15.0)
    # osim_template = OSIMSequence.a_pose(position = [0,0,1], name='osim_template')
    
    # c = (85/255, 85/255, 255/255, 1)

    # colors = vertex_colors_from_weights(weights=range(len(markers_pos_rel)), scale_to_range_1=True, alpha=1)[np.newaxis, :, :]
    # markers_rel_pc = PointClouds(markers_pos_rel, colors=colors, name='markers_rel')



    # #rajagopal unposed mesh
    # rajagopal_trimesh = trimesh.load('/home/kellerm/Data/OpenSim/rajagopal.obj', process=None)
    # rajagopal_mesh = Meshes(rajagopal_trimesh.vertices, rajagopal_trimesh.faces, name='Rajagopal')



    # markers_pos_rel, markers_pos_abs = get_markers_location(skel_osim)
    if marker_set == 'skin_set':
        # Reconstruct skin mesh form the markers 
        smpl_mesh = trimesh.load_mesh(cg.sample_smpl_mesh, process=False) # Load a SMPL mesh to get the faces
        
        vertices = osim_model_results.marker_trajectory[0]
        faces = smpl_mesh.faces
        rajagopal_skin_mesh = Meshes(vertices=vertices, faces=faces, position = [0,0,-1], name='Rajagopal_skin')
        v.scene.add(rajagopal_skin_mesh)

    # Display in the viewer

    # v.run_animations = True
    v.scene.camera.position = np.array([5.0, 2.5, 0.0])
    # v.scene.add(osim_template)
    # v.scene.add(markers_abs_pc, markers_rel_pc, rajagopal_mesh)

    return v

