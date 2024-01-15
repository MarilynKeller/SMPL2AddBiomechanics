import argparse
import json
import os
import pickle
import tqdm

from smpl2ab.markers.smpl2osim import Smpl2osim
from smpl2ab.markers.smpl_markers import SmplMarker
import smpl2ab.config as cg

from smpl2ab.measurements.measurements import BodyMeasurements
from smpl2ab.utils.smpl_utils import SMPL, load_smpl_seq, smpl_model_fwd

def create_subj_json_and_mesh(smpl_data):
    
    gender = smpl_data['gender']
    betas = smpl_data['betas']

    body_measurements = BodyMeasurements.from_smpl_params(gender, betas[:10])

    subject_json = {
        "subjectTags": ["healthy"], 
        "sex": str(gender), 
        "massKg": f"{body_measurements.compute_mass():0.2f}", 
        "heightM": f"{body_measurements.compute_height():0.2f}", 
        "email": "", 
        "skeletonPreset": "custom"}

    return subject_json


def create_data_folder(subject_name, subject_trials, output_folder, force_recompute=False, marker_set_name='bsm'):
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
    subject_json = create_subj_json_and_mesh(seq_data)
    print('Subject json: ', subject_json)
    json.dump(subject_json, open(subject_json_file, 'w'))
    print('Subject measurements saved as {}'.format(subject_json_file))
    
    smpl_model = SMPL(gender = seq_data['gender'])
    
    # Select the marker set to use
    if marker_set_name == "bsm":
        from markers.marker_sets import bsm_marker_set
        markers_dict = bsm_marker_set
        osso_segmentation = pickle.load(open(cg.bsm_osso_segmentation, 'rb'))
        rigging_method = 'gt'
    else:
        rigging_method = 'closest_rj_bone'
        raise ValueError(f"Unknown marker set: {marker_set_name}. To create a new marker set, add it as a dictionary in smpl_marker_dict.py")
    
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

    # markers_dict = smpl_manual_markers
    osim_model_path = cg.bsm_model_path
    assert os.path.exists(osim_model_path), f'Could not find {osim_model_path}'
    sos = Smpl2osim.from_files(markers_dict, osim_model_path, osso_segmentation, marker_set_name=marker_set_name, rigging_method=rigging_method)   
    # sos.generatse_osim(smpl_mesh_path=smpl_mesh_path, osso_mesh_path=osso_mesh_path, output_osim_path=output_osim_path, display=display)


    # Populate trial folder with synthetic mocap
    for seq in tqdm.tqdm(subject_trials):
        
        seq_name = seq.split('/')[-1].split('.')[0]
        
        # Create a folder for the sequence
        seq_folder = os.path.join(trial_folder, seq_name)
        os.makedirs(seq_folder, exist_ok=True)
        print(f'Created folder {seq_folder}')

        # Generate the synthetic markers for this sequence
        synth_mocap_file = os.path.join(seq_folder, 'markers.trc')
        if not os.path.exists(synth_mocap_file) or force_recompute:
           
            # Load the SMPL sequence as a dictionary
            seq_data = load_smpl_seq(seq)
            
            synthetic_markers = SmplMarker.from_smpl_data(smpl_data=seq_data, marker_set_name=marker_set_name, smpl_model=smpl_model)
            synthetic_markers.save_trc(synth_mocap_file)
            pickle.dump(synthetic_markers, open(synth_mocap_file, 'wb'))
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
    parser.add_argument('-o','--output_folder', type=str,help='Path to the output folder', default='./output')
    parser.add_argument('-f','--force_recompute', action='store_true', help='Force recomputing the synthetic markers')
    
    args = parser.parse_args()
    
    subject_trials = list_trials(args.smpl_seq_folder)
    subject_name = os.path.basename(args.smpl_seq_folder)
    
    print(f"Found {len(subject_trials)} trials for subject {subject_name}")
    
    create_data_folder(subject_name, subject_trials, args.output_folder, force_recompute=args.force_recompute)