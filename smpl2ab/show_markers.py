# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import argparse
import os
import numpy as np
import yaml

from aitviewer.renderables.osim import OSIMSequence
from aitviewer.renderables.smpl import SMPLSequence
from aitviewer.renderables.markers import Markers
from aitviewer.viewer import Viewer

import config as cg
from smpl2ab.markers.smpl_markers import SmplMarker

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--osim_path', type=str, default=cg.osim_model_path, help='Path to OpenSim model')
    parser.add_argument('--smpl_markers_path', type=str, default=cg.bsm_markers_on_smpl_path, help='Path to SMPL markers')
    
    args = parser.parse_args()
    
    to_display = []
    
    seq_smpl = SMPLSequence.t_pose()
    to_display.append(seq_smpl)
    
    markers_dict = yaml.load(open(args.smpl_markers_path, 'r'), Loader=yaml.FullLoader)
    synthetic_markers = SmplMarker(seq_smpl.vertices, markers_dict, fps=60, name='SMPL markers')
    markers_seq = Markers(synthetic_markers.marker_trajectory, markers_labels=synthetic_markers.marker_names, name='SMPL markers')
    to_display.append(markers_seq)
    

    osim_seq = OSIMSequence.a_pose(osim_path=args.osim_path ,color_skeleton_per_part=True, color_markers_per_part=True, position=[1,-0.5,0])   
    to_display.append(osim_seq)
    
    # Display in the viewer
    v = Viewer()
    v.scene.camera.position = np.array([0,0, 4.0])
    v.scene.add(*to_display)    
    
    
    v.run()
