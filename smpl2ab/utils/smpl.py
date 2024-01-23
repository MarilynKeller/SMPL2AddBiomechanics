# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 
 
import os
import pickle
import torch
import trimesh
import numpy as np
import smplx

from aitviewer.models.smpl import SMPLLayer
import config as cg

NUM_POSE_PARAMS = 72 # only use the first 72 joints (the rest are for the hands)

class SmplData():
    
    def __init__(self, poses, trans, betas, gender, mocap_framerate=None, use10betas=False, unpose_hands=False) -> None:
        if use10betas:
            betas = betas[:,:10]
            
        if unpose_hands:
            poses[:, -6:] = 0 #unpose hands
            
        self.poses = poses
        self.trans = trans
        self.betas = betas
        self._model = None
        self.mocap_framerate = mocap_framerate
        
        if isinstance(gender, list):
            assert len(set(gender)) == 1, "The data given to SMPL contains different genders, this is not support for now."
            gender = gender[0]
        self.gender = gender
            
    @property
    def model(self):
        if self._model is None:
            if len(self.betas.shape) == 1:
                num_betas = self.betas.shape[0]
            elif len(self.betas.shape) == 2:
                num_betas = self.betas.shape[1]
            else:
                raise ValueError(f"Betas should be a 1D or 2D array, got {self.betas.shape}")
            batch_size = self.poses.shape[0]
            self._model = smplx.create(cg.smpl_folder, model_type='smpl', gender=self.gender, num_betas=num_betas, batch_size=batch_size)
        return self._model

    @classmethod
    def from_amass(cls, amass_file:str, unpose_hands=False, use10betas=False, verbose=True):
        
        assert os.path.exists(amass_file), f'{amass_file} does not exist'
        
        if verbose:
            print(f"loading {amass_file}")
        
        if amass_file.endswith('.pkl'):
            smpl_data = pickle.load(open(amass_file, 'rb'))
            if 'YOGI' in amass_file:
                gender = 'female'
            else:
                gender = smpl_data['gender'] 
            trans = smpl_data['transl']
            poses = np.zeros((smpl_data['transl'].shape[0], 72))
            poses[:,:3] = smpl_data['global_orient']
            poses[:,3:72] = smpl_data['body_pose']
            betas = np.mean(smpl_data['betas'],axis=0)

        elif amass_file.endswith('.npz'):
            smpl_data = np.load(amass_file)

            gender = smpl_data['gender']
            if isinstance(gender, np.ndarray):
                gender = str(gender)

            poses = smpl_data['poses'][:, :NUM_POSE_PARAMS] 
            trans = smpl_data['trans']
            betas = smpl_data['betas']
            
            if use10betas:
                betas = betas[:,:10]
            gender = gender  
            
        else:
            raise ValueError(f"amass_file should be a .pkl or .npz file, got {amass_file}")
        
        if "mocap_framerate" in dict(smpl_data).keys():
            mocap_framerate = smpl_data['mocap_framerate']
        else:
            print(f'No mocap info in {amass_file}')
            mocap_framerate = 60 # default mocap_framerate
            
        if isinstance(mocap_framerate, np.ndarray):
            mocap_framerate = float(mocap_framerate)
        
        if unpose_hands:
            # unpose hands
            poses[:, -6:] = 0.0 # set the last 6 pose params to 0
        
        return cls( poses=poses, trans=trans, betas=betas, gender=gender, mocap_framerate=mocap_framerate)
    
    def get_vertices(self, is_amass=True):
      
        smpl_model = self.model
        
        if len(self.betas.shape) == 1:
            torch_betas = torch.FloatTensor(self.betas).unsqueeze(0).expand(self.poses.shape[0], -1)
        elif len(self.betas.shape) == 2:
            torch_betas = torch.FloatTensor(self.betas)
        else:
            raise ValueError(f"Betas should be a 1D or 2D array, got {self.betas.shape}")
        body_pose = torch.FloatTensor(self.poses[:,3:NUM_POSE_PARAMS])
        global_orient = torch.FloatTensor(self.poses[:, 0:3])
        transl = torch.FloatTensor(self.trans)
        
        smpl_output = smpl_model(betas=torch_betas, body_pose=body_pose, transl=transl, global_orient=global_orient)
        vertices = smpl_output.vertices
         
        if is_amass:
            # AMASS bodies are flipped, the joints were corrected before but the vertices need to be rotated to superimpose the joints
            vertices = vertices[:,:,[0,2,1]]
            vertices[:,:,[2]] =  - vertices[:,:,[2]]

        return vertices
    
    @property
    def faces(self):
        return self.model.faces_tensor
    

def amass2params(amass_file):
    """ from an amass file, load the model parameters
    @param amass_file: path to the amass file
    @return: smpl_poses, smpl_beta, gender
    """
    assert os.path.exists(amass_file), f'{amass_file} does not exist'
    print(f"loading {amass_file}")
    smpl_data = np.load(amass_file)
    smpl_poses = smpl_data['poses'][:, :NUM_POSE_PARAMS] 
    smpl_trans = smpl_data['trans']
    smpl_poses = np.hstack([smpl_trans, smpl_poses])
    smpl_beta = smpl_data['betas']
    
    if isinstance(gender, np.ndarray):
        gender = str(gender)
    
    return smpl_poses, smpl_beta, gender


def load_mean_smpl(device='cpu'):

    smpl_layer = SMPLLayer(model_type='smpl', gender='male', device=device)
    poses = torch.zeros([1, smpl_layer.bm.NUM_BODY_JOINTS * 3], device=device)
    betas = torch.zeros([1, smpl_layer.num_betas], device=device)

    verts, joints = smpl_layer.fk(poses, betas)

    return trimesh.Trimesh(verts.detach().cpu().numpy()[0], smpl_layer.bm.faces, process=False)

