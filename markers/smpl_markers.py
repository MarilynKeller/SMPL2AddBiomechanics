import numpy as np
from smpl_utils import smpl_model_fwd
import nimblephysics as nimble

class SmplMarker():
    """ Class that, given a SMPL sequence, generates a marker sequence from specific SMPL vertices."""
        
    def __init__(self, verts, markers_dict, fps, name):
        self.verts = verts
        self.set = markers_dict  
        self.fps = fps
        self.name = name

        markers_smpl_indices = list(markers_dict.values())

        self.marker_trajectory = verts[:, markers_smpl_indices]
        self.marker_names = list(markers_dict.keys())


    def save_trc(self, path):

        markerTimesteps = []
        for frame_id in range(len(self.marker_trajectory)):
            markers_dict = {}
            for marker_id, marker_name in enumerate(self.marker_names):
                markers_dict[marker_name] = self.marker_trajectory[frame_id][marker_id]
            markerTimesteps.append(markers_dict)

        timestamps = list(np.arange(0, self.marker_trajectory.shape[0]/self.fps, 1/self.fps))
        timestamps = list( 1/self.fps * np.arange(0, len(markerTimesteps)))

        assert(len(timestamps) == len(markerTimesteps))
        nimble.biomechanics.OpenSimParser.saveTRC(path, timestamps, markerTimesteps)
        
    @classmethod
    def from_smpl_data(cls, smpl_data, marker_set, smpl_model, fps=None):
        """ Create a marker sequence from a SMPL sequence"""
        
        # Select the marker set to use
        if marker_set == "bsm":
            from markers.marker_sets import bsm_marker_set
            markers_dict = bsm_marker_set
        else:
            raise ValueError(f"Unknown marker set: {marker_set}. To create a new marker set, add it as a dictionary in smpl_marker_dict.py")

        # FPS for the sequence
        if fps is None and "fps" in smpl_data.keys():
            fps = smpl_data['fps']
        elif fps is not None:
            # Use the value provided as argument
            pass
        else:
            raise ValueError("The motion sequence FPS was not found in smpl_data (expected a key 'fps'), please provide it as argument")
        
        # Per frame SMPL vertices to generate the markers from
        verts = smpl_model_fwd(smpl_model, smpl_data)

        return cls(verts, markers_dict, fps, name=marker_set)
    