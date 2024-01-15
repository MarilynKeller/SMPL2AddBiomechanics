import os
import pickle

import numpy as np
import torch
import smplx
import config as cg

def load_smpl_seq(smpl_seq_path, include_verts=False):

    if not os.path.exists(smpl_seq_path):
        raise Exception('Path does not exist: {}'.format(smpl_seq_path))
    
    if smpl_seq_path.endswith('.pkl'):
        data_dict = pickle.load(open(smpl_seq_path, 'rb'))
    
    elif smpl_seq_path.endswith('.npz'):
        data_dict = np.load(smpl_seq_path)
        data_dict = {key: data_dict[key] for key in data_dict.keys()} # convert to python dict
        
        # In some npz, the gender type happens to be: array('male', dtype='<U4'). So we convert it to string
        if not isinstance(data_dict['gender'], str):
            data_dict['gender'] = str(data_dict['gender'])
            
        if data_dict['poses'].shape[1] == 156:
            # Those are SMPL+H poses, we remove the hand poses to keep only the body poses
            poses = np.zeros((data_dict['poses'].shape[0], 72))
            poses[:, :72-2*3] = data_dict['poses'][:, :72-2*3] # We leave params for SMPL joints 22 and 23 to zero 
            data_dict['poses'] = poses
        
    for key in ['trans', 'poses', 'betas', 'gender']:
        assert key in data_dict.keys(), f'Could not find {key} in {smpl_seq_path}. Available keys: {data_dict.keys()})'
        
    out_dict = {}
    out_dict['trans'] = data_dict['trans']
    out_dict['poses'] = data_dict['poses']
    out_dict['betas'] = data_dict['betas']
    out_dict['gender'] = data_dict['gender']
    
    if "mocap_framerate" in data_dict.keys():
        out_dict['fps'] = data_dict['mocap_framerate']
        
    return out_dict


def smpl_model_fwd(smpl_model, smpl_data, device='cpu'):
    # Run a SMPL forward pass to get the SMPL numpy body vertices from  a numpy SMPL sequence dictionary
    to_torch = lambda x: torch.from_numpy(x).float().to(device)
    
    poses_smpl = to_torch(smpl_data['poses'])
    trans_smpl = to_torch(smpl_data['trans'])
    betas_smpl = to_torch(smpl_data['betas'][:smpl_model.num_betas]).expand(trans_smpl.shape[0], -1)
    
    # Run a SMPL forward pass to get the SMPL body vertices
    smpl_output = smpl_model(betas=betas_smpl, body_pose=poses_smpl[:,3:], transl=trans_smpl, global_orient=poses_smpl[:,:3])
    verts = smpl_output.vertices
    return verts
   
    
def SMPL(gender, num_betas=10): 
    smpl= smplx.create(
            model_path=cg.smpl_folder,
            gender=gender,
            num_betas=num_betas,
            model_type='smpl')
    return smpl
    