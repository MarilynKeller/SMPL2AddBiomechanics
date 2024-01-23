# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 
 
import argparse
from typing import Optional
import os

from aitviewer.renderables.markers import Markers
from aitviewer.renderables.osim import OSIMSequence
from aitviewer.utils.vtp_to_ply import convert_meshes
from aitviewer.viewer import Viewer

import smpl2ab.config as cg

def display_model_in_viewer(osim: Optional[str] = None, mot: Optional[str] = None, fps: int = 5,
                            color_parts: bool = False, color_markers: bool = False, mocap: Optional[str] = None,
                            joints: bool = False) -> None:
    """Load an OpenSim model and a motion file and display it in the viewer.

    Args:
        osim (str, optional): Path to the osim file. If no file is specified, the default Rajagopal gait model will be loaded from nimble. Defaults to None.
        mot (str, optional): Path to the motion file. Defaults to None.
        fps (int, optional): Generating the meshes for all the frames can take a long time and a lot of memory. Use a low fps to avoid this problem. Defaults to 30.
        color_parts (bool, optional): Color the skeleton by parts, as defined in the .osim file. Defaults to False.
        color_markers (bool, optional): Each marker is attached to a skeleton part. This option colors the markers to match the parent part color. Defaults to False.
        mocap (str, optional): If a Mocap file is specified, display the markers mocap with their labels. For now, only the .c3d format is supported. Defaults to None.
        joints (bool, optional): Show model joints as spheres. Defaults to False.
    """
    
    
    if osim is not None:
        # Check that a folder named Geometry exists in the same folder as the osim file.
        osim_dir = os.path.dirname(osim)
        geom_dir = os.path.join(osim_dir, "Geometry")
        if not os.path.exists(geom_dir):
            print("Geometry folder does not exist in the same folder as the osim file. Please add a Geometry folder containing the skeleton .vtp, .obj or .ply file in the same folder as the provided .osim. ")
            exit(1)
            
        # Check that the Geometry folder contains at least a file of type .ply
        ply_files = [f for f in os.listdir(geom_dir) if f.endswith(".ply")]
        if len(ply_files) == 0:
            print("Geometry folder does not contain any .ply files.")
            print("The provided folder meshes will be converted to .ply. Press a key to continue or CTRL-C to abort.")

            # Convert the provided meshes to .ply
            convert_meshes(geom_dir, geom_dir)

    if mot is None:
        osim_seq = OSIMSequence.a_pose(osim, name='OpenSim template', 
                                        show_joint_angles=joints, 
                                        color_skeleton_per_part=color_parts, 
                                        color_markers_per_part=color_markers,
                                        )
        
    else:
        osim_seq = OSIMSequence.from_files(osim, mot, 
                                        show_joint_angles=joints, 
                                        color_skeleton_per_part=color_parts, 
                                        color_markers_per_part=color_markers,
                                        fps_out=fps, position = [-1,0,0], ignore_fps=True)
    
    
    v = Viewer()
    v.scene.add(osim_seq)
    
    if mocap is not None:
        # check that the mocap file is in .c3d format
        assert mocap.endswith(".c3d"), "Mocap file must be in .c3d format."
        marker_seq = Markers.from_c3d(mocap, fps_out=fps, color=[0,255,0,255])
        v.scene.add(marker_seq)
        
    v.playback_fps = fps

    v.run()
    

class KinHelper():
    """ Load a Open Sim model to list the bone groups and markers to help buildng the correspondances with osim"""
    def __init__(self, osim_model_path, mot=0): 
        self.osim_file = osim_model_path 
        if mot is not None:
            self.mot = mot
        else:
            self.mot = cg.osim_sample_motion
        
    def display_osim(self):
        display_model_in_viewer(osim = self.osim_file, mot = None,
                            color_parts = False, color_markers = False,
                            joints = True)

    def display_osim_motion(self):
        
        display_model_in_viewer(osim = self.osim_file, mot = self.mot,
                            color_parts = True, color_markers = True,
                            joints = True)  
              
    def print_joint_labels_dict(self):
        """Example KinOsim.print_pose_labels('/ps/project/rib_cage_breathing/TML/Data/AddBiomechanics/Models/fused_shoulder_spine/fused_shoulder_spine.osim')"""
        import nimblephysics as nimble
        osim = nimble.biomechanics.OpenSimParser.parseOsim(self.osim_file)
        node_names = [n.getName() for n in osim.skeleton.getBodyNodes()]
        nodes_dict = {i:n for i,n in enumerate(node_names)}
        print(nodes_dict)
         
    def print_params_labels_dict(self):
        """Example KinOsim.print_pose_labels('/ps/project/rib_cage_breathing/TML/Data/AddBiomechanics/Models/fused_shoulder_spine/fused_shoulder_spine.osim')"""
        import nimblephysics as nimble
        osim = nimble.biomechanics.OpenSimParser.parseOsim(self.osim_file)
        node_names = [q.getName() for q in osim.skeleton.getDofs()]
        params_label_dict = {i:n for i,n in enumerate(node_names)}
        print(params_label_dict)
        
    def print_marker_labels_list(self):
        """Example KinOsim.print_marker_labels('/ps/project/rib_cage_breathing/TML/Data/AddBiomechanics/Models/fused_shoulder_spine/fused_shoulder_spine.osim')"""
        import nimblephysics as nimble
        osim = nimble.biomechanics.OpenSimParser.parseOsim(self.osim_file)
        markers_labels = [ml for ml in osim.markersMap.keys()]
        markers_labels.sort()
        print(markers_labels)
        
    def print_marker_rigging_list(self):
        """Example KinOsim.print_marker_labels('/ps/project/rib_cage_breathing/TML/Data/AddBiomechanics/Models/fused_shoulder_spine/fused_shoulder_spine.osim')"""
        import nimblephysics as nimble
        osim = nimble.biomechanics.OpenSimParser.parseOsim(self.osim_file)
        markers_labels = [ml for ml in osim.markersMap.keys()]
        marker_parents = []
        for marker in markers_labels:
            parent = osim.markersMap[marker][0].getName()
            marker_parents.append(parent)
            print(marker, parent)
        # markers_labels.sort()
      
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser('Analyse a .osim file to build corresp with SMPL. By default show the default .osim in use in the project', 
                                     usage='python kin_models/kin_helpers.py -o /ps/project/rib_cage_breathing/TML/Data/AddBiomechanics/Models/fused_shoulder_spine/fused_shoulder_spine.osim')
    parser.add_argument('--osim_file', '-o', type=str, help='Path to the .osim file', default=cg.osim_model_path)
    parser.add_argument('--display', '-D', action='store_true', help='Display the osim model')
    parser.add_argument('--motion', '-M', action='store_true', help='Display the osim model doing a motion')
    parser.add_argument('--mot', type=str, default=None, help='Motion file to viszualize')
    
    parser.add_argument('--joints', '-j', action='store_true', help='Print joints name of the model')
    parser.add_argument('--params', '-p', action='store_true', help='Print list of parameters of the model')
    parser.add_argument('--markers', '-m', action='store_true', help='Print markers of the model')
    parser.add_argument('--markers_rigging', '-mr', action='store_true', help='Print markers rigging of the model: (marker name, parent bone))')
    
    args = parser.parse_args()
    
    kin_helper = KinHelper(args.osim_file, args.mot)
    
    if args.display:
        kin_helper.display_osim()
        
    if args.motion or args.mot is not None:
        kin_helper.display_osim_motion()
    
    if args.joints:
        print('OSIM model joints:')
        kin_helper.print_joint_labels_dict()  
    
    if args.markers:
        print('OSIM model markers:')  
        kin_helper.print_marker_labels_list()
    
    if args.params:
        print('OSIM model parameters:')  
        kin_helper.print_params_labels_dict()
        
    if args.markers_rigging:
        print('OSIM model markers rigging:')  
        kin_helper.print_marker_rigging_list()