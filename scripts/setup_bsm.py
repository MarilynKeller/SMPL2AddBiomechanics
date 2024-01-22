# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import argparse
import os
import shutil
import zipfile


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Place the bsm model at the proper place')
    
    parser.add_argument('skel_zip', type=str, default='/path/to/skel_models_v1.x.zip')
    
    args = parser.parse_args()
    
    model_folder = 'models'
    skel_zip = args.skel_zip
    skel_zip_name = os.path.basename(skel_zip)
    skel_zip_name = skel_zip_name.replace('.zip', '')
    
    # Unzip the BSM model
    print('Unzipping BSM model...')
    with zipfile.ZipFile(skel_zip, 'r') as zip_ref:
        zip_ref.extractall(model_folder)

    os.makedirs(os.path.join(model_folder, 'bsm'), exist_ok=True)
    
    # Rename the SMPL files
    print('Reorganizing BSM files...')
    shutil.move(os.path.join(model_folder, f'{skel_zip_name}/bsm.osim'),
                os.path.join(model_folder, 'bsm', 'bsm.osim'))
    
    if not os.path.exists(os.path.join(model_folder, 'bsm', 'Geometry')):
        shutil.move(os.path.join(model_folder, f'{skel_zip_name}/Geometry'),
                os.path.join(model_folder, 'bsm', 'Geometry'))
    
    motion_dir = os.path.join(model_folder, 'bsm', 'sample_motion')
    os.makedirs(motion_dir, exist_ok=True)
    if not os.path.exists(os.path.join(motion_dir, '01')):
        shutil.move(os.path.join(model_folder, f'{skel_zip_name}/sample_motion/'),
            os.path.join(motion_dir, '01'))
    
    skel_seq = os.path.join(model_folder, 'bsm', 'sample_motion/01/01_01_poses_skel.pkl')
    if os.path.exists(skel_seq):
        os.remove(skel_seq)

    print("Cleaning up...")
    shutil.rmtree(os.path.join(model_folder, f'{skel_zip_name}'))
    print("Done!")
    print(" The BSM models was placed in models/.")