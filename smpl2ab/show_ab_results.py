# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import argparse
import os
import numpy as np
import yaml

from aitviewer.renderables.osim import OSIMSequence
from aitviewer.renderables.smpl import SMPLSequence
from aitviewer.renderables.markers import Markers
from aitviewer.viewer import Viewer
from aitviewer.models.smpl import SMPLLayer

import config as cg
from smpl2ab.markers.smpl_markers import SmplMarker
from smpl2ab.utils.smpl_utils import load_smpl_seq

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--osim_path', type=str, help='Path to OpenSim model (.osim)')
    parser.add_argument('--mot_path', type=str,  help='Path to OpenSim motion (.mot)')
    parser.add_argument('--smpl_motion_path', type=str,  help='Path to SMPL motion')
    parser.add_argument('--smpl_markers_path', type=str, default=cg.bsm_markers_on_smpl_path, help='Path to SMPL markers')
    parser.add_argument('--body_model', help='Body model to use (smpl or smplx)', default='smpl', choices=['smpl', 'smplx'])
    parser.add_argument('--gender', type=str, default=None, choices=['female', 'male', 'neutral'])
    parser.add_argument('--z_up', action='store_true', help='Set the z axis up')
     
    args = parser.parse_args()
    
    to_display = []
    
    if args.body_model == 'smpl':
        body_model = 'smplh'
    else :
        body_model = args.body_model
        
    if args.gender is None:
        smpl_data = load_smpl_seq(args.smpl_motion_path)
        gender = smpl_data['gender']
    else:
        gender = args.gender
        
    smpl_layer = SMPLLayer(model_type=body_model, gender=gender)
    args = parser.parse_args()
    
    to_display = []
    
    fps = 30 

    # Load SMPL motion
    seq_smpl = SMPLSequence.from_amass(
        smpl_layer = smpl_layer,
        npz_data_path=os.path.join(args.smpl_motion_path), # AMASS/CMU/01/01_01_poses.npz
        fps_out=fps,
        name=f"{args.body_model.upper()} motion",
        show_joint_angles=False,
        z_up=args.z_up
    )
    to_display.append(seq_smpl)
   
    # Load result OpenSim motion
    osim_seq = OSIMSequence.from_files(osim_path=args.osim_path, 
                                       mot_file=args.mot_path, 
                                       name=f'OpenSim skeleton motion', 
                                       fps_out = fps,
                                       color_skeleton_per_part=False, 
                                       show_joint_angles=False, 
                                       is_rigged=False,
                                       ignore_geometry=True,
                                       z_up=args.z_up)
    
    to_display.append(osim_seq)


    # Load SMPL markers
    markers_dict = yaml.load(open(args.smpl_markers_path, 'r'), Loader=yaml.FullLoader)
    synthetic_markers = SmplMarker(seq_smpl.vertices, markers_dict, fps=fps, name='Markers')
    markers_seq = Markers(synthetic_markers.marker_trajectory, markers_labels=synthetic_markers.marker_names, 
                          name='SMPL markers',
                          color=(0, 1, 0, 1),
                          z_up=args.z_up)
    to_display.append(markers_seq)
    

    # Display in the viewer
    v = Viewer()
    v.run_animations = True
    v.scene.camera.position = np.array([10.0, 2.5, 0.0])
    v.scene.add(*to_display)
    
    if seq_smpl is not None:
        v.lock_to_node(seq_smpl, (2, 0.7, 2), smooth_sigma=5.0)
    v.playback_fps = fps
    
    v.run()