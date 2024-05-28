# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import argparse
import logging
import os
import numpy as np
import trimesh
import pickle as pkl

from aitviewer.renderables.meshes import Meshes
from aitviewer.viewer import Viewer
from aitviewer.models.smpl import SMPLLayer
from aitviewer.renderables.markers import Markers
from aitviewer.utils.colors import vertex_colors_from_weights

from smpl2ab.markers.mapping import apply_mapping, apply_mapping_barycentric, compute_mapping, compute_mapping_barycentric
from smpl2ab.utils.osso import get_osso_submesh_util 
from smpl2ab.utils.smpl import load_mean_smpl
from smpl2ab.utils.osim_aug import OsimAug
from smpl2ab.utils.smpl2bsm import smpl2osim_corresp

import config as cg

def flatten(l):
    return [item for sublist in l for item in sublist]

class SSMarkerTransfer():

    osso_rj_segmentation = cg.bsm_osso_segmentation 
    output_folder = cg.smpl_markers_folder

    def __init__(self, markers_p, osso_p, osso_rj, markers_label, skin_indices, mapping_method, rigging_method):

        # meshes
        self.markers_p = markers_p #smpl posed 
        self.osso_p = osso_p # osso posed
        self.osso_rj = osso_rj# osso dismembered to be aligned with Rajagopal unposed

        # markers relative to osso
        self.markers_osso_d = None
        self.markers_osso_idx = None
        self.m_osso_vect = None

        self.m_label = markers_label 
        self.skin_indices = skin_indices # index of the skin vertices each marker is attached to

        # markers_rigging
        self.markers_rigging = None

        # markers relative to unposed rajagopal
        self.markers_rj = None

        self.mapping_method = mapping_method
        self.rigging_method = rigging_method

        self.osso_rj_seg = pkl.load(open(cg.bsm_osso_segmentation, 'rb'))
        
        self._femur_demo_cache = {'bone': 'femur_l'} # Store the femur mesh and associated markers


    def compute_rel_markers(self, method=None):

        if method is None:
            method = self.mapping_method
        assert method in ['relative', 'absolute', 'per_part'] # baricentric should not be used 

        if method == 'absolute':
            # Get the offset vector in world coordinates
            self.markers_osso_d, self.markers_osso_idx  = trimesh.proximity.ProximityQuery(self.osso_p).vertex(self.markers_p) #vertex index on osso
            self.m_osso_vect = self.markers_p - self.osso_p.vertices[self.markers_osso_idx] #skeleton to marker vectors

            # Compute markers location on unposed RJ
            self.markers_rj = self.reconstruct_markers(self.osso_rj, from_osim=True)
            # Check that all the markers got reconstructed
            assert self.markers_rj.shape == self.markers_p.shape

        elif method == 'barycenter':
            self.mesh_tri_index, self.baricentric_verts_coords = compute_mapping_barycentric(self.markers_p, self.osso_p) 
 
            # TEST
            # m_rec = apply_mapping_barycentric(self.baricentric_verts_coords, self.mesh_tri_index, self.osso_p)
            # assert(np.linalg.norm(m_rec - self.markers_p))
            self.markers_rj = apply_mapping_barycentric(self.baricentric_verts_coords, self.mesh_tri_index, self.osso_rj)
            self.markers_rj[:, [0,2]] = self.markers_rj[:, [2,0]]
            self.markers_rj[:, 2] = -self.markers_rj[:, 2]
            print('WARNING: flipping the result')

            # return self.markers_osso_idx, self.markers_triangle_id, self.m_osso_vect_rel
            
        elif method == 'relative':
            self.mesh_tri_index, self.m_vect_rel = compute_mapping(self.markers_p, self.osso_p) 

            # Compute markers location on unposed RJ
            self.markers_rj = self.reconstruct_markers(self.osso_rj)
            # Check that all the markers got reconstructed
            assert self.markers_rj.shape == self.markers_p.shape

        elif method == 'per_part':
            self.mesh_tri_index = np.zeros((self.markers_p.shape[0]), dtype=np.int32)
            self.m_vect_rel = np.zeros((self.markers_p.shape[0], 3)) 

            #self.markers_rigging[m_label] : rajagopal_bone_name
            for bone_name in self.osso_rj_seg.keys():
              
                # Failed markers 
                failed_markers = [ml for ml in self.m_label if ml not in self.markers_rigging.keys()]
                print(f'WARNING: The following markers were not rigged: {failed_markers}, please add those makers manually to the osim model')
                # sys.exit()
                
                # Find indices of the markers attached to the bone
                valid_markers = [ml for ml in self.m_label if (ml not in failed_markers)]
                markers_mask = [self.m_label.index(ml) for ml in valid_markers if (self.markers_rigging[ml] == bone_name)]
                
                markers = self.markers_p[markers_mask]

                if len(markers) == 0:
                    print('WARNING: No markers attached to bone {}'.format(bone_name))
                    continue

                submesh, submesh_global_faces_index = self.get_osso_submesh(bone_name)

                # Compute mapping to go from this bone triangles to the markers  
                submesh_tri_index, m_vect_rel = compute_mapping(markers, submesh) 

                # Convert the submesh triangle index to the corresponfing full mesh triangle index
                mesh_tri_index = submesh_global_faces_index[submesh_tri_index]
            
                # Save the mapping for this bone markers
                self.mesh_tri_index[markers_mask] = mesh_tri_index
                self.m_vect_rel[markers_mask] = m_vect_rel
                
                if bone_name == self._femur_demo_cache['bone']:
                    self._femur_demo_cache['mesh'] = self.get_osso_unposed_submesh(bone_name)[0]
        
            self.markers_rj = self.reconstruct_markers(self.osso_rj, from_osim=False)

        else:
            raise ValueError('Unknown mapping method')


    def reconstruct_markers(self, osso_mesh, from_osim=True, method=None):
        """
        Reconstruct the markers from the osso mesh. If you want to reconstruct the markers from rajagopal, set from_osim to True.
        """
        if method is None:
            method = self.mapping_method
        if method == 'absolute':
            # From a mesh with osso topology, no matter its pose, reconstruct the markers location
            m_reconstructed = osso_mesh.vertices[self.markers_osso_idx ]  \
                            + self.m_osso_vect
                            # + m_osso_dist[:,np.newaxis] * osso_rajagopal_trimesh.vertex_normals[m_osso_idx] 

        elif method == 'barycenter':
            m_reconstructed = apply_mapping_barycentric(self.baricentric_verts_coords, self.mesh_tri_index, osso_mesh)

        elif method in ['relative', 'per_part']:
            m_reconstructed = apply_mapping(self.mesh_tri_index, self.m_vect_rel, osso_mesh)
        
        else:
            raise ValueError('Unknown mapping method')

        if from_osim:
            m_reconstructed[:, [0,2]] = m_reconstructed[:, [2,0]]
            m_reconstructed[:, 2] = -m_reconstructed[:, 2]
            print('WARNING: flipping the result')
                
        return m_reconstructed


    def rig_markers_to_osim_bones_using_closest_point(self):
        """
        For each skin marker, pick the closest osso vertex. loook in a precomputed dict 
        to which rajagopal bone this vertex corresponds, then rig the marker to this bone
        """
        
        # Load osso segmentation into RJ bones, this segmentation was created with blender
        vertex_rigging_dict =  pkl.load(open(SSMarkerTransfer.osso_rj_segmentation, 'rb'))
        # logging.info([(k, v) for (k,v) in vertex_rigging_dict.items()])
        self.markers_osso_d, self.markers_osso_idx  = trimesh.proximity.ProximityQuery(self.osso_p).vertex(self.markers_p) #vertex index on osso

        #For each marker and its corresponding vertex on osso
        self.markers_rigging = {}
        for mi, (m_osso_index, m_label) in enumerate(zip(self.markers_osso_idx, self.m_label)):
            m_osso_part = None
            for bone_name, bone_vertices in vertex_rigging_dict.items():
                if m_osso_index in bone_vertices:
                    m_osso_part = bone_name
                    break
            if m_osso_part is None:
                print(f'WARNING: Could not find the bone corresponding to marker {m_label}')
            else:
                self.markers_rigging[m_label] = vertex_rigging_dict[m_osso_part]
        logging.info([(k, v) for (k,v) in self.markers_rigging.items()])


    def rig_markers_to_osim_bones_using_smpl_rig(self):
        """
        For each skin marker, pick the bone of smpl that has the highest skinning weight, 
        look in a dict to find the corresponfing rajagopal bone, and rig the marker to this rajagopal bone.
        """
        # Load dict mapping smpl bone index to rajagopal bone name
        from auto_markers.rajagopal_smpl_rigging import smpl_to_raj_dict 
        
        smpl_layer = SMPLLayer(model_type='smpl', gender='neutral', device='cpu')
        skining_weights = smpl_layer.bm.lbs_weights.numpy()
    
        #For each marker 
        self.markers_rigging = {}

        #Build a dict of distance queries
        distance_queries = {}

        ambiguous_bone = flatten([bones for bones in smpl_to_raj_dict.values() if len(bones)>1])
        for bone_name in ambiguous_bone:
            submesh, _ = self.get_osso_submesh(bone_name)
            distance_queries[bone_name] = trimesh.proximity.ProximityQuery(submesh)
        
        # iterate through (marker_coords, index in smpl, marker label)
        # m_smpl_pos : 3D coordinates of the marker
        smpl2osim_list = [v for v in smpl2osim_corresp.values()] # smpl2osim_list[0] contains ['pelvis'], this is the list of osim nodes corresponding to smpl joint 0
        for mi, (m_smpl_pos, m_smpl_index, m_label) in enumerate(zip(self.markers_p, self.skin_indices, self.m_label)):
            smpl_bone_index = np.argmax(skining_weights[m_smpl_index]) # Find to which smpl bone the marker is attached
            rajagopal_bone = self.find_rigging(vertex=m_smpl_pos, bone_candidates = smpl2osim_list[smpl_bone_index], bone_candidates_queries=distance_queries) # Get the corresponding rajagopal bone
            self.markers_rigging[m_label] = rajagopal_bone # Rig the marker to the rajagopal bone
        [print('\n',(k, v)) for (k,v) in self.markers_rigging.items()]
        
        
    def rig_markers_to_osim_bones_using_gt(self):
        '''Use the rigging specified in the reference osim model. Here the rigging is predifined and not computed from the mesh'''
        osim_aug = OsimAug()
        marker_labels, rajagopal_bones = osim_aug.list_markers_rigging()
        self.markers_rigging = {}
        self._femur_demo_cache['marker_indices'] = []
        for mi, (markers_label, rajagopal_bone) in enumerate(zip(marker_labels, rajagopal_bones)):
            self.markers_rigging[markers_label] = rajagopal_bone
            
            if rajagopal_bone == self._femur_demo_cache['bone']:
                self._femur_demo_cache['marker_indices'].append(mi)
            
        assert [k for k in self.markers_rigging.keys()] == sorted(self.markers_rigging.keys()), f'The created markers dict is not sorted properly'
    
    def find_rigging(self, vertex, bone_candidates, bone_candidates_queries):
        """Given a vertex and candidate bones, find the bone that is closest to the vertex"""
        assert len(bone_candidates) > 0, 'The mapping dict from smpl to rajagopal should at least list one bone'
        if len(bone_candidates) == 1:
            return bone_candidates[0]
        else:
            # Among the bone candidates, find the one that is closest to the vertex
            closest_bone = bone_candidates[0]
            closest_bone_dist = +np.inf
            for bone_name in bone_candidates:
                bone_proximity_query = bone_candidates_queries[bone_name]
                _, bone_dist, _ = bone_proximity_query.on_surface(vertex[np.newaxis,:])
                # print(f"{bone_name} distance: {bone_dist[0]}")
                if bone_dist[0] < closest_bone_dist:
                    closest_bone_dist = bone_dist
                    closest_bone = bone_name

            return closest_bone


    def get_osso_submesh(self, bone_name):
        return get_osso_submesh_util(self.osso_p, self.osso_rj_seg, bone_name)

    def get_osso_unposed_submesh(self, bone_name):
        return get_osso_submesh_util(self.osso_rj, self.osso_rj_seg, bone_name)

    def rig_markers_to_osim_bones(self):
        method = self.rigging_method
        if method == 'closest_rj_bone':
            self.rig_markers_to_osim_bones_using_closest_point()
        elif method == 'heavy_smpl_bone':
            self.rig_markers_to_osim_bones_using_smpl_rig()
        elif method == 'gt':
            self.rig_markers_to_osim_bones_using_gt()
        else:
            raise ValueError(f'Unknown rigging method {method}, should be one of closest_rj_bone, heavy_smpl_bone, gt')


    def export_correspondances(self, marker_set_name):
        
        # Rig each marker wrt rajagopal bones
        self.rig_markers_to_osim_bones()

        #Compute mapping from osso to rajagopal unposed
        self.compute_rel_markers()

        assert (self.markers_rigging is not None) and (self.markers_rj is not None), 'Data not created yet, run self.rig_markers_to_osim_bones()'

        # Save correspondances between markers and osso bones
        out_folder = os.path.join(SSMarkerTransfer.output_folder, marker_set_name)
        os.makedirs(out_folder, exist_ok=True)
        loc_file = out_folder + f'/{marker_set_name}_locations.npy'
        rig_file = out_folder + f'/{marker_set_name}_rigging.npy'
        np.save(loc_file, self.markers_rj)
        pkl.dump(self.markers_rigging, open(rig_file, 'wb'))

        print(f'The markers correspondance files were created as : \n\t{loc_file}\n\t{rig_file}')

    @classmethod
    def load_rigging(cls, marker_set_name):
        """Return location, rigging_dict of the marker set"""
        out_folder = os.path.join(SSMarkerTransfer.output_folder, marker_set_name)
        loc_file = out_folder + f'/{marker_set_name}_locations.npy'
        rig_file = out_folder + f'/{marker_set_name}_rigging.npy'
        location = np.load(loc_file)
        rigging_dict = pkl.load(open(rig_file, 'rb'))
        print(f'Loading correspondance files : \n\t{loc_file}\n\t{rig_file}')
        return location, rigging_dict


    def visualize(self, smpl_trimesh):
        
        skin_mesh_color = (179/255, 120/255, 179/255, 0.5)

        reconstructed_p = self.reconstruct_markers(self.osso_p, from_osim=False)
        
        # Set offset to display the different meshes
        dx = np.array([1, 0.0, 0.0])
        
        # Load unposed rajagopal mesh for visualization
        rj = trimesh.load_mesh(cg.osim_unposed, process=False)
        
        # To show only some markers
        if False:
            M = 500
            self.markers_p = self.markers_p[:M]
            reconstructed_p = reconstructed_p[:M]
            self.markers_rj = self.markers_rj[:M]
        
        # Set the color of the markers
        marker_index_colors = vertex_colors_from_weights(weights=range(len(self.markers_p)), scale_to_range_1=True, alpha=1)#[np.newaxis, :, :]

        disp_offset = 0*dx
        smpl_mesh = Meshes(smpl_trimesh.vertices, smpl_trimesh.faces, name='SMPL mesh', position=disp_offset, color = skin_mesh_color)
        osso_mesh = Meshes(self.osso_p.vertices, self.osso_p.faces, name='OSSO mesh', position=disp_offset)
        
        # if len(self.markers_p) > 200:
        #     Points = PointClouds
        # else:
        Points = Markers
        # markers_labels= [f'{i}' for i in range(len(self.markers_p))]
        markers_labels = self.m_label
         
        markers_pc = Points(self.markers_p[None,...], markers_labels=markers_labels,  position=disp_offset, colors=marker_index_colors, name='Markers on SMPL')
        m_osso_rec_pc = Points(reconstructed_p[None,...], markers_labels=markers_labels, position=disp_offset, colors=marker_index_colors, name='Markers reconstructed from OSSO')
              
        # Demo bone
        print(self._femur_demo_cache['marker_indices'])
        demo_bone_mesh = self._femur_demo_cache['mesh']
        demo_bone_marker_indices = self._femur_demo_cache['marker_indices']
        
        disp_offset = 1*dx
        markers_labels_demo = [f'{i}' for i in demo_bone_marker_indices]
        m_osso_pc_demo = Points(self.markers_rj[demo_bone_marker_indices][None,...], markers_labels=markers_labels_demo, position=disp_offset, colors=marker_index_colors[demo_bone_marker_indices], name='Markers reconstructed from osso_rajagopal bone')
        osso_mesh_demo = Meshes(demo_bone_mesh.vertices, demo_bone_mesh.faces, name='Demo bone', position=disp_offset)

        # Markers anchor on OSSO
        disp_offset = 2*dx
        self.markers_osso_d, self.markers_osso_idx  = trimesh.proximity.ProximityQuery(self.osso_p).vertex(self.markers_p)
        m_osso_pc = Points(self.osso_p.vertices[self.markers_osso_idx][None,...], markers_labels=markers_labels, position=disp_offset, colors=marker_index_colors, name='Markers anchor on OSSO')
        osso_mesh2 = Meshes(self.osso_p.vertices, self.osso_p.faces, name='OSSO mesh 2', position=disp_offset)
        # line_renderable = LinesWithGeometryShader.from_start_end_points(self.markers_rj[np.newaxis,:], markers_p, r_base=0.005, name='lines')
        # link_skels = LinesWithGeometryShader.from_start_end_points(osso_p.vertices[np.newaxis,:], osso_rj.vertices,r_base=0.002, name='lines')

        disp_offset = 3*dx
        rj_mesh = Meshes(rj.vertices, rj.faces, name='Rajagopal unposed mesh', position=disp_offset)
        m_rajagopal_pc2 = Points(self.markers_rj[None,...], markers_labels=markers_labels, position=disp_offset, colors=marker_index_colors, name='Markers reconstructed from osso_rajagopal')

        disp_offset = 4*dx
        # m_rajagopal_pc = PointClouds(self.markers_rj, position=p, color=[1, 0, 0, 1], name='markers reconstructed from osso_rajagopal')
        rj_mesh1 = Meshes(rj.vertices, rj.faces, name='Rajagopal mesh', position=disp_offset, color=(0.7,0.7,1,1))
        osso_rj_mesh = Meshes(self.osso_rj.vertices, self.osso_rj.faces, name='OSSO Rajagopal', position=disp_offset, color=(0.5,0.5,0,1))

        # Display in viewer.
        v = Viewer()

        v.scene.add(
                    smpl_mesh, 
                    osso_mesh, 
                    m_osso_rec_pc, 
                    markers_pc,
                    
                    m_osso_pc, 
                    osso_mesh2,
                    
                    # m_rajagopal_pc, 
                    osso_rj_mesh,
                    rj_mesh1,
                    # line_renderable,
                    # link_skels

                    m_rajagopal_pc2,
                    rj_mesh,
                    
                    m_osso_pc_demo,
                    osso_mesh_demo,
                    )
        v.run()



def get_vertices(smpl_trimesh, marker_set_name):

    if marker_set_name == 'skin_set':
        labels = [f'{i:04d}' for i in range(len(smpl_trimesh.vertices))]
        skin_indices = [i for i in range(len(smpl_trimesh.vertices))]
        return smpl_trimesh.vertices, labels, skin_indices

    elif marker_set_name == 'smpl_cmu':
        smpl_cmu_dict = pkl.load(open(cg.smpl_cmu_file, 'rb'))
        skin_indices = np.hstack(smpl_cmu_dict.values())
        markers_label = [k for k in smpl_cmu_dict.keys()]
        return smpl_trimesh.vertices[skin_indices], markers_label, skin_indices
    else:
        raise ValueError(f'Unknown marker set name {marker_set_name}')




if __name__ == '__main__':

    parser = argparse.ArgumentParser('Create a marker set relative to rajagopal. Given vertices defined on SMPL, \
        this script will store the location and rigging of each marker wrt rajagopal')
    parser.add_argument('-s','--marker_set_name', type=str, default='skin_set', help='Name of the marker set', choices=['skin_set', 'smpl_cmu'])
    parser.add_argument('-D','--display', action='store_true', default=False, help='Display the result in a viewer')
    parser.add_argument('-rm','--rigging_method', type=str, required=True, help='Rigging method', choices=['heavy_smpl_bone', 'closest_rj_bone'])
    parser.add_argument('-mm', '--mapping_method', type=str, required=True, help='Mapping method', choices=['relative', 'absolute', 'barycenter', 'per_part'])

    args = parser.parse_args()
    marker_set_name = args.marker_set_name

    pose = 'lying_down'
    # pose = 'tpose'

    if pose == 'tpose': 
        # Create markers from SMPL vertices
        smpl_trimesh = load_mean_smpl()
        markers_p, markers_label, skin_indices = get_vertices(smpl_trimesh, marker_set_name)

        # Load the corresponding OSSO skeleton
        osso_p = trimesh.load_mesh('/home/kellerm/Data/OSSO/template.obj', process=False)
        osso_rj = trimesh.load_mesh(cg.osso_rj_unposed, process=False)

    if pose == 'lying_down':
        smpl_trimesh = trimesh.load_mesh('/home/kellerm/Data/AMASS_OSSO/OSSO_raw/CMU/01/01_03_poses/template/star_lying.ply', process=False)
        markers_p, markers_label, skin_indices = get_vertices(smpl_trimesh, marker_set_name)

        # Load the corresponding OSSO skeleton
        osso_p = trimesh.load_mesh('/home/kellerm/Data/AMASS_OSSO/OSSO_raw/CMU/01/01_03_poses/template/skel_lying.ply', process=False)
        # osso_p = to_template_topo(osso_p)
        osso_rj = trimesh.load_mesh(cg.osso_rj_unposed, process=False)

    mt = SSMarkerTransfer(markers_p, osso_p, osso_rj, markers_label, skin_indices, rigging_method=args.rigging_method, mapping_method=args.mapping_method)

    mt.export_correspondances(marker_set_name)

    if args.display:
        mt.visualize(smpl_trimesh)



