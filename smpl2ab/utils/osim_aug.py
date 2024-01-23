# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import os
import trimesh
import nimblephysics as nimble
import numpy as np

from aitviewer.renderables.osim import OSIMSequence
from aitviewer.viewer import Viewer

import config as cg


class OsimAug():
    
    def __init__(self, osim_path=cg.osim_model_path) -> None:
        
        assert os.path.exists(osim_path), f'{osim_path} does not exist'
        self.osim_path = osim_path
        self.osim = nimble.biomechanics.OpenSimParser.parseOsim(osim_path)
        self.skeleton = self.osim.skeleton
        
        self.node_names = [n.getName() for n in self.osim.skeleton.getBodyNodes()]
        self._meshes_dict = {}
        self._indices_dict = {}
        self._generate_meshes_dict() 
           
           
    def load_motion(self, motion_path):
        "Load the motion and return the pose parameters as a numpy array"
        mot: nimble.biomechanics.OpenSimMot = nimble.biomechanics.OpenSimParser.loadMot(
                    self.osim.skeleton, motion_path)

        motion = np.array(mot.poses.T)  
        return motion  

 
 
    def create_template(self):
        part_meshes = []
        for node_name in self.node_names:
            mesh = self.meshes_dict[node_name]
            # assert mesh, "No mesh for node: {}".format(node_name)
            if mesh is None:
                print( "WARNING: No mesh for node: {}".format(node_name))
            if mesh:
                part_meshes.append(mesh)
        self.template = trimesh.util.concatenate(part_meshes)
        
    
    def get_bone_submesh(self, node:str) -> trimesh.Trimesh:
        """ Return the submesh for the given node """
        return self.meshes_dict[node]
    
    
    
    def get_markers_location(self):
        ############ Get markers locations
        markers_abs = self.skeleton.getMarkerMapWorldPositions(self.osim.markersMap)
        print(f"Number of markers: {len(markers_abs)}")

        markers_pos_rel = np.vstack([m[1] for m in self.osim.markersMap.values()])
        markers_pos_abs = np.vstack(markers_abs.values())

        ################# Get markers bone attachments
        # body_nodes = [m[0] for m in osim.markersMap.values()]
        # for bn in body_nodes:
        #     print(bn.getName())

        return markers_pos_rel, markers_pos_abs
    
    
    
    def list_bone_infos(self):
        """ List the bones and the corresponding meshes"""
        node_names = [n.getName() for n in self.osim.skeleton.getBodyNodes()]
        print('Number of nodes:', len(node_names))
        print(node_names)

        # List the meshes constituting each part
        for node_name in node_names:
            print(node_name)
            for shape_nodes in self.osim.skeleton.getBodyNode(node_name).getShapeNodes():
                print( '\t', shape_nodes.getShape().getMeshPath().split('/')[-1])
                
    def list_markers_rigging(self):
        ''' Return two lists: the first one contains the markers labels, the second one the parent bone of each marker'''
        markers_labels = [ml for ml in self.osim.markersMap.keys()]
        marker_parents = []
        for marker in markers_labels:
            parent = self.osim.markersMap[marker][0].getName()
            marker_parents.append(parent)
            print(marker, parent)
        return markers_labels, marker_parents 
        
        
    def _generate_meshes_dict(self):
        """ Output a dictionary giving for each bone, the attached mesh"""

        current_index = 0
        self.indices_dict = {}
        self.meshes_dict = {}

        node_names = self.node_names
        for node_name in node_names:
            mesh_list = []
            body_node = self.osim.skeleton.getBodyNode(node_name)
            # print(f' Loading meshes for node: {node_name}')
            num_shape_nodes = body_node.getNumShapeNodes()
            if num_shape_nodes == 0:
                print(f'WARNING:\tNo shape nodes listed for  {node_name}')
            for shape_node_i in range(num_shape_nodes):
                shape_node = body_node.getShapeNode(shape_node_i)
                submesh_path = shape_node.getShape().getMeshPath()
                # Get the scaling for this meshes
                scale = shape_node.getShape().getScale()
                offset = shape_node.getRelativeTranslation()
                # Load the mesh
                try:
                    submesh = trimesh.load_mesh(submesh_path, process=False)
                    # print(f'Loaded mesh {submesh_path}')
                except Exception as e:
                    print(e)
                    print(f'WARNING:\tCould not load mesh {submesh_path}')
                    submesh = None
                    continue
                
                if submesh is not None:
                    trimesh.repair.fix_normals(submesh)
                    trimesh.repair.fix_inversion(submesh)
                    trimesh.repair.fix_winding(submesh)

                    # Scale the bone to match .osim subject scaling
                    submesh.vertices[:] = submesh.vertices * scale
                    submesh.vertices[:] += offset
                    # print(f'submesh_path: {submesh_path}, Nb vertices: {submesh.vertices.shape[0]}')
                    mesh_list.append(submesh)

            # Concatenate meshes
            if mesh_list:
                node_mesh = trimesh.util.concatenate(mesh_list)
                self.indices_dict[node_name] = (current_index, current_index + node_mesh.vertices.shape[0])
                current_index += node_mesh.vertices.shape[0]
            else:
                node_mesh = None
                print("\t WARNING: No submesh for node:", node_name)
                self.indices_dict[node_name] = (current_index, current_index )
            
            # Add to the dictionary
            self.meshes_dict[node_name] = node_mesh
        print(self.meshes_dict)
        
    def get_joints_location(self):
        node_names = self.node_names
        joints_loc = np.zeros((len(node_names), 3))
        for ni, node_name in enumerate(node_names):

            transfo = self.osim.skeleton.getBodyNode(node_name).getTransform()
            joints_loc[ni] = transfo.translation()
        return joints_loc
    
    def get_bone_scales(self):
        node_names = self.node_names
        scales = np.zeros((len(node_names), 3))
        for ni, node_name in enumerate(node_names):
            body_node = self.osim.skeleton.getBodyNode(node_name)
            for shape_node in body_node.getShapeNodes():
                scales[ni] = shape_node.getShape().getScale()
        return scales
                
        
        
    def visualize(self):     
    
        v = Viewer()

        osim_model_results = OSIMSequence.a_pose(osim_path=self.osim_path, name='result_osim', color_markers_per_part=True, color_skeleton_per_part=True)
        v.scene.add(osim_model_results)
        
        v.run()
              
    def get_joints(self):
        node_names = self.node_names
        joints_loc = np.zeros((len(node_names), 3))
        for ni, node_name in enumerate(node_names):
            transfo = self.osim.skeleton.getBodyNode(node_name).getTransform()
            joints_loc[ni] = transfo.translation()
        return joints_loc
    
    def find_rot_axis(self):
        node_names = self.node_names
        joints_loc = np.zeros((len(node_names), 3))
        for ni, node_name in enumerate(node_names):
            transfo = self.osim.skeleton.getBodyNode(node_name).getTransform()
            joints_loc[ni] = transfo.translation()
        return joints_loc
        
    def get_joints_flip_axis_map(self):
        joints = self.osim.skeleton.getJoints()
        for j in joints:
            try:
                print(j.getFlipAxisMap())
            except:
                print('No coordinate')
   
        
if __name__ == '__main__':
    
    oa = OsimAug()  
    oa.get_joints_flip_axis_map()