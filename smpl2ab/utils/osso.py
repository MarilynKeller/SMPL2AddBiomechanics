# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import os
from typing import Iterable, Union
import numpy as np
import trimesh
import pickle as pkl
import config as cg
import nimblephysics as nimble

from smpl2ab.markers.mapping import get_submesh

""" Class to handle the segmentation of an osso mesh in osim nodes"""
class OssoAug(trimesh.Trimesh):
    
    def __init__(self, osso_mesh_path=None) -> None:
        
        osim_file = cg.osim_model_path
        osso_segmentation_file = cg.bsm_osso_segmentation

        if osso_mesh_path == None:
            osso_mesh_path = cg.osso_inf_template
        assert os.path.exists(osso_mesh_path), f'{osso_mesh_path} does not exist'
                
        self._osim = nimble.biomechanics.OpenSimParser.parseOsim(osim_file)
        self._osso_segmentation = pkl.load(open(osso_segmentation_file, 'rb'))

        self.mesh = trimesh.load(osso_mesh_path, process=False, maintain_order=True)
        self.bone_names = list(self._osso_segmentation.keys())
    
   
    def get_bone_indices(self, bone_names: list) -> list:
        # Return a list of vertices indices for the given bone names
        return get_bone_indices_utils(self._osso_segmentation, bone_names)
 

    def get_osso_submesh(self, bone_name: list) -> Iterable[Union[trimesh.Trimesh, np.array]]:
        """Given a mesh with osso template topology and a rajagopal bone name, return the submesh of osso corresponding to this bone"""
        if bone_name not in self.bone_names:
            raise ValueError('Bone name not found')
        return get_osso_submesh_util(self.mesh, self._osso_segmentation, bone_name)
  
  
def get_bone_indices_utils(osso_rj_seg, bone_names: list) -> list:
    # Return a list of vertices indices for the given bone names
    if isinstance(bone_names, list):
        bone_indices = []
        for bone in bone_names:
            bone_indices.append(np.array(osso_rj_seg[bone]))
        bone_indices = np.concatenate(bone_indices)
    else:
        bone_indices = np.array(osso_rj_seg[bone_names])
    return bone_indices  
    
# This function needs to be independant of the class because I use it in building .osim from unposed osso 
def get_osso_submesh_util(osso_mesh, osso_rj_seg, bone_name):     
    """Given a mesh with osso template topology and a rajagopal bone name, return the submesh of osso corresponding to this bone"""

    # Build a binary mask of the osso vertices to keep
    bone_indices = get_bone_indices_utils(osso_rj_seg, bone_name)

    # Build a submesh of the bone
    new_verts, new_faces, bool_faces, vertex_ids = get_submesh(osso_mesh.vertices, osso_mesh.faces, verts_retained=bone_indices)
    
    # Compute offset between the original mesh vertices and the created submesh. Should be zero.
    diff = np.max(np.abs(osso_mesh.vertices[bone_indices] - new_verts))
    if diff != 0 : 
        print(f'diff : {diff}, there is a mismatch between the original vertices and the submesh')
 
    # Check output
    print('bone_name: ', bone_name, '\tvertices_nb: ', np.count_nonzero(bone_indices), '\tnew_verts_nb: ', new_verts.shape[0], '\tnew_faces_nb: ', new_faces.shape[0], '\tfaces_bool: ', np.count_nonzero(bool_faces))
    
    submesh = trimesh.Trimesh(vertices=new_verts, faces=new_faces, process=False, maintain_order=True)
    submesh_global_faces_index = np.where(bool_faces)[0] # For each face of the submesh, gices its index in the original mesh
    
    # Check that i can recover the original faces/triangle indices
    assert np.all((submesh.triangles - osso_mesh.triangles[bool_faces]) == 0)

    if new_verts.shape[0] != len(bone_indices):
        raise ValueError('The submesh for {bone_name} does not have the expected number of verts. {new_verts.shape[0]} vs {len(bone_indices)} expected')

    return submesh, submesh_global_faces_index
