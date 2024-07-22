# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

from typing import NewType
import torch
import trimesh
import yaml
import numpy as np
import smplx
import config as cg

Tensor = NewType('Tensor', torch.Tensor)

class BodyMeasurements:
    """ Class that, given betas and gender of the SMPL model, computes the height and weight of a smpl body."""       
    
    DENSITY = 985
    
    def __init__(self, verts, faces, triangles, model_type='smpl'):

        self.verts =  verts
        self.faces = faces
        self.triangles = triangles
        self.model_type = model_type
        
        self.print_measurements()


    @classmethod
    def from_smpl_params(cls, gender, betas, model_type='smpl', device='cpu'):
        smpl= smplx.create(
            model_path=cg.smpl_folder,
            gender=gender,
            num_betas=betas.shape[0],
            model_type=model_type
        ).to(device)

        # T pose
        pose = torch.FloatTensor(np.zeros(smpl().body_pose.shape[1]), device=device).unsqueeze(0)
        betas = torch.FloatTensor(betas, device=device).unsqueeze(0)

        # Apply T pose and betas
        body = smpl(poses_body = pose, betas = betas[:])
        triangles = body.vertices[:, smpl.faces_tensor]
        return cls(body.vertices, smpl.faces_tensor, triangles, model_type) 
    
    @classmethod
    def from_mesh_file(cls, mesh_path):
        mesh = trimesh.load(mesh_path, process=False)
        vertices = torch.tensor(mesh.vertices, dtype=torch.float32).unsqueeze(0)
        faces = torch.tensor(mesh.faces, dtype=torch.int64)
        triangles = vertices[:, faces]
        return cls(vertices, faces, triangles)

    @property
    def body_mesh(self):
        mesh = trimesh.Trimesh(self.verts.detach().numpy()[0], self.faces)
        return mesh

    def compute_height(self):
        ''' Compute the height using the heel and the top of the head
        Code adapted from Lea Muller (lea.muller@tuebingen.mpg.de), https://github.com/muelea/shapy/tree/master/measurements, 
        '''
        
        mesh_vertices_path = cg.mesh_vertices_smplx_path if self.model_type=='smplx' else cg.mesh_vertices_path
        with open(mesh_vertices_path, 'r') as f:
            meas_vertices = yaml.safe_load(f)

        head_top = meas_vertices['HeadTop']
        left_heel = meas_vertices['HeelLeft']

        left_heel_bc = left_heel['bc']
        left_heel_face_idx = left_heel['face_idx']

        left_heel_bc = torch.tensor(left_heel['bc'], dtype=torch.float32)
        head_top_bc = torch.tensor(head_top['bc'], dtype=torch.float32)

        head_top_face_idx = head_top['face_idx']

        head_top_tri = self.triangles[:, head_top_face_idx]
        head_top = (head_top_tri[:, 0, :] * head_top_bc[0] +
                    head_top_tri[:, 1, :] * head_top_bc[1] +
                    head_top_tri[:, 2, :] * head_top_bc[2])
        head_top = (
            head_top_tri * head_top_bc.reshape(1, 3, 1)
        ).sum(dim=1)
        left_heel_tri = self.triangles[:, left_heel_face_idx]
        left_heel = (
            left_heel_tri * left_heel_bc.reshape(1, 3, 1)
        ).sum(dim=1)

        return torch.abs(head_top[:, 1] - left_heel[:, 1]).item()


    def compute_mass(self):
        ''' Computes the mass from volume and average body density
        Code adapted from Lea Muller (lea.muller@tuebingen.mpg.de), https://github.com/muelea/shapy/tree/master/measurements, 
        '''
        x = self.triangles[:, :, :, 0]
        y = self.triangles[:, :, :, 1]
        z = self.triangles[:, :, :, 2]
        volume = (
            -x[:, :, 2] * y[:, :, 1] * z[:, :, 0] +
            x[:, :, 1] * y[:, :, 2] * z[:, :, 0] +
            x[:, :, 2] * y[:, :, 0] * z[:, :, 1] -
            x[:, :, 0] * y[:, :, 2] * z[:, :, 1] -
            x[:, :, 1] * y[:, :, 0] * z[:, :, 2] +
            x[:, :, 0] * y[:, :, 1] * z[:, :, 2]
        ).sum(dim=1).abs() / 6.0
        mass = volume * self.DENSITY
        return mass.item()


    def print_measurements(self):
        ''' Prints the measurements
        '''
        print(f"Height: {self.compute_height()}")
        print(f"Weight: {self.compute_mass()}")


