# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import os
import pickle
import shutil
import numpy as np
import trimesh
import nimblephysics as nimble

from aitviewer.renderables.meshes import Meshes
from aitviewer.renderables.osim import OSIMSequence
from aitviewer.viewer import Viewer

from smpl2ab.markers.marker_transfer import SSMarkerTransfer
import config as cg

class OsimEditor:
    
    def __init__(self, input_osim_path, marker_set_name) -> None:
        
        self.input_osim_path = input_osim_path
        self.marker_set_name = marker_set_name
        
        self.osim = nimble.biomechanics.OpenSimParser.parseOsim(self.input_osim_path)
        
        
    def export_osim(self, output_osim_path):
        
        # Create is_anatomical dictionary
        is_anatomical = {}
        anatomical_list = self.osim.anatomicalMarkers
        for marker_id, marker_obj in self.osim.markersMap.items():
            if marker_id in anatomical_list:
                is_anatomical[marker_id] = True
                print(f"Marker {marker_id} marked as anatomical")
            else:
                is_anatomical[marker_id] = False
            
        marker_dict = create_osim_marker_dict(self.osim, self.marker_set_name)
        

        output_osim_path_gen = output_osim_path.replace('.osim', '_gen.osim')
        os.makedirs(os.path.dirname(output_osim_path_gen), exist_ok=True)
        nimble.biomechanics.OpenSimParser.replaceOsimMarkers(self.input_osim_path, marker_dict, is_anatomical, output_osim_path_gen)
        print(f'Osim model with marker set {self.marker_set_name} saved as ', output_osim_path)

        shutil.move(output_osim_path_gen, output_osim_path)

        # Save an array of which marker is rigged to each bone for faster visualisation in aitviewer
        rigging_path = output_osim_path.replace('.osim', f'_rigging.pkl')
        save_rigging(self.osim, marker_dict, rigging_path)
        print(f"Rigging saved to {rigging_path}")
        
        # Copy geometry from source to target osim
        source_geom = '/'.join(self.input_osim_path.split('/')[:-1]) + '/Geometry'
        target_geom = '/'.join(output_osim_path.split('/')[:-1]) + '/Geometry'

        # Copy the geometry folder to cg.current_model_cluster_folder
        # import ipdb; ipdb.set_trace()
        # shutil.copytree(source_geom, cg.osim_geometry_folder, dirs_exist_ok=True, symlinks=True) # Copy the whole geometry folder
        
        # if os.path.islink(target_geom):
        #     os.unlink(target_geom)
        # os.symlink(cg.osim_geometry_folder, target_geom) # Create a symbolic link of the geometry folder
        # print(f"Geometry copied from {source_geom} to {target_geom}")
        
       
       

def create_osim_marker_dict(osim, marker_set_name):
    """ To modify the osim markers, we need to create a dictionary of the form: markers_dict =  'M0': ('pelvis', [0, 1, 2]).
    This function does it """

    locations, rigging = SSMarkerTransfer.load_rigging(marker_set_name)
     
    osim.markersMap
    markersMap = {} # This is the marker obj used by the osim obj
    # To save the markers as a .obj file, I need to provide it under the format  markers_dict =  'M0': ('pelvis', [0, 1, 2])}
    marker_dict = {} 
    skeleton_body_nodes = osim.skeleton.getBodyNodes()
    skeleton_body_nodes_names = [n.getName() for n in skeleton_body_nodes]
    for mi, (m_label, node_name_list) in enumerate(rigging.items()):
        # if vertex_id>10:
        #     break
        # Ignore markers that are not rigged to any bone
        # if not node_name_list:
        #     continue
        assert(node_name_list), f"{m_label} is not rigged to any bone. In Blender, make sure to rig each vertex to a bone."
        if isinstance(node_name_list, list) :
            node_name = node_name_list[0]
        else:
            node_name = node_name_list
        location = locations[mi]
        try:
            bodyNode = skeleton_body_nodes[skeleton_body_nodes_names.index(node_name)]
        except Exception as e:
            print(e)
        
        # For some anatomical markers, we placed them very  carefully ao we keep tho osim onces
        osim_location = np.array(osim.markersMap[m_label][1])
        
        # For those anatomical markers, we trust better the manual definition than OSSO
        if m_label in ['LAKI', 'LELB', 'LFIN', 'LHEE', 'LHEB', 'LKNE', 'LTOE', 
                       'LFPI', 'LWRA', 'LWRB', 'RANK', 'RSCA', 'RELB', 'RFIN', 
                       'RHEE', 'RHEB', 'RTOE', 'RFPI', 'RWRA', 'RWRB', 'RHPI', 'RHME', 
                       'RHTH', 'LHPI', 'LHME', 'LHTH', 'LKNI', 'RKNI', 
                       'LANK', 'RAKI', 'RTOS', 'RTOP', 'LTOS', 'LTOP', 'LSCA', 'LELS', 'LELSO', 'RELS', 'RELSO',
                       'LTIA', 'LTIB', 'LTIC', 'RTIA', 'RTIB', 'RTIC',
                       'RFRM', 'LFRM',
                       'C7', 'LUMB', 'LUMC', 'CLAV',
                       ]:
            print(f"Correcting anatomical marker {m_label} from {location} to {osim_location}. Difference: {np.linalg.norm(location-osim_location)}")
            
            location = osim_location
                
        # If osso places the markers inside the body, use OSIM default markers. This corrects  for petruding sternum and hips
        if m_label in ['LFWT', 'RFWT', 'RIBR','RIMR', 'STRN', 'STRM', 'LIBR','LIMR']:
            if np.linalg.norm(osim_location)>np.linalg.norm(location):
                location = osim_location

        # For the .osim dictionary
        marker_dict[str(m_label)]  = (node_name, location)
        markersMap[str(m_label)] =  (bodyNode, location)

    #Replace existing markers by skin markers
    osim.markersMap = markersMap

    return marker_dict


def save_rigging(osim, marker_dict, rigging_path):
    ''' Save the rigging as a pickle file.
        File format: { 'bones': ['pelvis', 'femur_l', ...]
                       'per_marker_rigging': [0, 0, 0, 1, 1, 2, 1, ... ] } # Index of the bone each marker is rigged to
    '''
    # We have marker_dict[str(m_label)]  : (node_name, location)
    bone_names = [n.getName() for n in osim.skeleton.getBodyNodes()]
    rigging_dict = {
        'bones': bone_names,
        'per_marker_rigging': [bone_names.index(k[0]) for k in marker_dict.values()]
    }

    with open(rigging_path, 'wb') as f:
        pickle.dump(rigging_dict, f)
    



def visualize_osim(osim_path, marker_set):


    v = Viewer()
    osim_model_results = OSIMSequence.a_pose(osim_path=osim_path, name='result_osim', color_markers_per_part=True, color_markers_per_index=False, color_skeleton_per_part=True)
    v.scene.add(osim_model_results)

    # markers_pos_rel, markers_pos_abs = get_markers_location(skel_osim)
    if marker_set == 'skin_set':
        # Reconstruct skin mesh form the markers 
        smpl_mesh = trimesh.load_mesh(cg.sample_smpl_mesh, process=False) # Load a SMPL mesh to get the faces
        
        vertices = osim_model_results.marker_trajectory[0]
        faces = smpl_mesh.faces
        rajagopal_skin_mesh = Meshes(vertices=vertices, faces=faces, position = [0,0,-1], name='Rajagopal_skin')
        v.scene.add(rajagopal_skin_mesh)

    # Display in the viewer
    v.scene.camera.position = np.array([5.0, 2.5, 0.0])

    return v

    
def corect_osim(src_path, dst_path):
    import xml.etree.ElementTree as ET

    # parse the XML file
    tree = ET.parse(src_path)
    root = tree.getroot()

    # Loop through all the Marker elements
    for marker in root.findall(".//Marker"):
        # Get the body element
        body = marker.find("body")
        # Create a new SocketParentFrame element
        socket_parent_frame = ET.Element("socket_parent_frame")
        # Set its text to the desired value
        socket_parent_frame.text = "/bodyset/" + body.text
        # Replace the body element with the new SocketParentFrame element
        # marker.replace(body, socket_parent_frame)
        marker.remove(body)
        # Insert the socket_parent_frame element in its place
        marker.insert(0, socket_parent_frame)


    # write modified XML to file
    tree.write(dst_path)
    print(f'Corrected .osim saved as {dst_path}')
    

