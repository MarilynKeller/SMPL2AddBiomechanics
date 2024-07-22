# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 
 
import argparse
import json
import os
import pickle
import tqdm
import yaml

from smpl2ab.markers.smpl2osim import Smpl2osim
from smpl2ab.markers.smpl_markers import SmplMarker
from smpl2ab.measurements.measurements import BodyMeasurements
from smpl2ab.utils.smpl_utils import SMPL, load_smpl_seq, smpl_model_fwd
import smpl2ab.config as cg

def create_subj_json_and_mesh(smpl_data, model_type='smpl'):
    
    gender = smpl_data['gender']
    betas = smpl_data['betas']

    body_measurements = BodyMeasurements.from_smpl_params(gender, betas[:10], model_type=model_type)

    subject_json = {
        "subjectTags": ["healthy"], 
        "sex": str(gender), 
        "massKg": f"{body_measurements.compute_mass():0.2f}", 
        "heightM": f"{body_measurements.compute_height():0.2f}", 
        "email": "", 
        "skeletonPreset": "custom"}

    return subject_json


def create_data_folder(subject_name, subject_trials, output_folder, osim_model_path=None, marker_dict_path=None, body_model='smpl', use_osso=False, force_recompute=False, display=True, no_confirm=False):
    """ Create the AdBiomechanics input data folder for a SMPL subject.
    Those data can be feed into AddBiomechanics to align an OpenSim Skeleton model to it.
    @param subject_name: name of the subject
    @param subject_trials: list of the paths of the SMPL sequences for this subject
    @param output_folder: path to output the generated data
    """

    # Create subject folder
    subject_folder = os.path.join(output_folder, subject_name)
    os.makedirs(subject_folder, exist_ok=True)

    # Create trial folder
    trial_folder = os.path.join(subject_folder, 'trials')
    if not os.path.exists(trial_folder):
        os.makedirs(trial_folder, exist_ok=True)

    # Generate subject_json with body measurements
    subject_json_file = os.path.join(subject_folder, '_subject.json')
    seq_data = load_smpl_seq(subject_trials[0]) # Load the first sequence
    subject_json = create_subj_json_and_mesh(seq_data, model_type=body_model)
    print('Subject json: ', subject_json)
    json.dump(subject_json, open(subject_json_file, 'w'))
    print('Subject measurements saved as {}'.format(subject_json_file))
    
    smpl_model = SMPL(gender = seq_data['gender'], model_type=body_model)
    
    # Select the marker set to use
    if osim_model_path == cg.osim_model_path and (marker_dict_path == cg.bsm_markers_on_smpl_path or marker_dict_path == cg.bsm_markers_on_smplx_path):
        # Use the default BSM model and marker set
        # from markers.marker_sets import bsm_marker_set
        # marker_set_name='bsm_smpl' if not use_smplx else 'bsm_smplx'
        if body_model == 'smpl':
            marker_set_name = 'bsm_smpl'
        elif body_model == 'smplx':
            marker_set_name = 'bsm_smplx'
        else:
            raise ValueError(f"Unknown body model {body_model}")
        
        with open(marker_dict_path, 'r') as yaml_file:
            markers_dict = yaml.safe_load(yaml_file)
        osso_segmentation = pickle.load(open(cg.bsm_osso_segmentation, 'rb'))
        rigging_method = 'gt'
    else :
        # Use a custom OpenSim model and marker set
        assert marker_dict_path is not None, 'If a custom osim model is used, a custom marker set must be provided'
        with open(marker_dict_path, 'r') as yaml_file:
            markers_dict = yaml.safe_load(yaml_file)
        print(f"WARNING: You are using a custom OpenSim model. If your model has different bone meshes than BSM, the marker transfer will fail. We advice to not use the --osso option in this case.")
        if not no_confirm:
            input('Press Ctrl+C to abort. If you wish to proceed anyway, press any other key ...')
        osso_segmentation = pickle.load(open(cg.bsm_osso_segmentation, 'rb'))
        marker_set_name = 'custom'
        rigging_method = 'closest_rj_bone'
        print(f'Using custom marker set: {markers_dict}')
        
  
    if use_osso:
        # Generate an BSM model with the corresponding markers
        osso_folder = os.path.join(subject_folder, 'osso')
        osso_mesh_path = os.path.join(osso_folder, 'skel_lying.ply')
        smpl_mesh_path = os.path.join(osso_folder, 'star_lying.ply') # smpl and star have the same mesh topology so they can be used interchangeably
        if not os.path.exists(osso_mesh_path):

            # Add vertices and face to seq_data
            verts = smpl_model_fwd(smpl_model, seq_data)
            faces = smpl_model.faces_tensor
            seq_data.update({'verts': verts, 'faces': faces})
            
            from osso.utils.fit_osso import fit_osso
            print('Fitting OSSO model to the SMPL mesh...')
            # This fit will be used to deduce the personalized offset between the bones and the skin markers
            fit_osso(seq_data, osso_folder, display=True)
            assert os.path.exists(osso_mesh_path), f'Could not find {osso_mesh_path}'
            assert os.path.exists(smpl_mesh_path), f'Could not find {smpl_mesh_path}'

        assert os.path.exists(osim_model_path), f'Could not find {osim_model_path}'
        sos = Smpl2osim.from_files(markers_dict, osim_model_path, osso_segmentation, marker_set_name=marker_set_name, rigging_method=rigging_method)   
        output_osim_path = os.path.join(subject_folder, 'bsm.osim')
        sos.generate_osim(smpl_mesh_path=smpl_mesh_path, osso_mesh_path=osso_mesh_path, output_osim_path=output_osim_path, display=display)


    # Populate trial folder with synthetic mocap
    for seq in tqdm.tqdm(subject_trials):
        
        seq_name = seq.split('/')[-1].split('.')[0]
        
        # Generate the synthetic markers for this sequence
        synth_mocap_file = os.path.join(trial_folder, seq_name+'.trc')

        if not os.path.exists(synth_mocap_file) or force_recompute:
           
            # Load the SMPL sequence as a dictionary
            seq_data = load_smpl_seq(seq)
            
            synthetic_markers = SmplMarker.from_smpl_data(smpl_data=seq_data, marker_set_name=marker_set_name, markers_dict=markers_dict, smpl_model=smpl_model)
            synthetic_markers.save_trc(synth_mocap_file)
            del synthetic_markers

            print(f'Generated synthetic markers as {synth_mocap_file}.')
        else:
            print(f'Synthetic markers already exist at {synth_mocap_file}. Not recomputing.')
            
  

def list_trials(smpl_seq_folder):
    """ Given a folder containing SMPL sequences, list the sequences paths """
    trials = []
    for trial in os.listdir(smpl_seq_folder):
        if trial.endswith('.npz') or trial.endswith('.pkl'):
            trials.append(os.path.join(smpl_seq_folder, trial))
    return trials


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Given a SMPL sequence, generate the input data for AddBiomechanics')
    
    parser.add_argument('-i', '--smpl_seq_folder', type=str, help='Path to the SMPL sequences folder', 
                        default='/is/cluster/work/mkeller2/Data/TML/AMASS/CMU/11')
    parser.add_argument('--osim', type=str, help='Path to a OpenSim model (.osim file)', default=cg.osim_model_path)
    parser.add_argument('--marker_dict', type=str, help='Path to a marker dictionary (.yaml file)', default=None)
    parser.add_argument('-o','--output_folder', type=str,help='Path to the output folder', default='./output')
    parser.add_argument('-F','--force_recompute', action='store_true', help='Force recomputing the synthetic markers')
    parser.add_argument('--body_model', help='Body model to use (smpl or smplx). If the sequences you use are SMPLH, use SMPL. If they are SMPLX, use SMPLX', default='smpl', choices=['smpl', 'smplx'])
    parser.add_argument('--osso', action='store_true', help='Gerenate personalized markers on the BSM template with the proper offsets to the bones, using OSSO.')
    parser.add_argument('-D','--display', action='store_true', help='If OSSO is used, display the result of the marker transfer.')
    parser.add_argument('--no_confirm', action='store_true', help="Do not ask for confirmation after warning.")
    
    args = parser.parse_args()
    
    subject_trials = list_trials(args.smpl_seq_folder)
    subject_name = os.path.basename(args.smpl_seq_folder)
    
    if len(subject_trials) == 0:
        raise ValueError(f"No SMPL sequences found at {args.smpl_seq_folder}")
    print(f"Found {len(subject_trials)} trials for subject {subject_name}")

    marker_dict_path = args.marker_dict
    # marker_dict_path = [marker_dict_path if marker_dict_path else [cg.bsm_markers_on_smplx_path if args.smplx else cg.bsm_markers_on_smpl_path][0]][0]
    if  marker_dict_path is None:
        # marker_dict_path = cg.bsm_markers_on_smplx_path if args.smplx else cg.bsm_markers_on_smpl_path
        if args.body_model == 'smplx':
            marker_dict_path = cg.bsm_markers_on_smplx_path
        elif args.body_model == 'smpl':
            marker_dict_path = cg.bsm_markers_on_smpl_path
        else:
            raise ValueError(f"Unknown body model {args.body_model}")

    create_data_folder(subject_name=subject_name, 
                       subject_trials=subject_trials, 
                       output_folder=args.output_folder, 
                       osim_model_path=args.osim, 
                       marker_dict_path=marker_dict_path, 
                       body_model=args.body_model, 
                       use_osso=args.osso, 
                       force_recompute=args.force_recompute, 
                       display=args.display,
                       no_confirm=args.no_confirm)