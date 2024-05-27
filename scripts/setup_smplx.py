# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import argparse
import os
import shutil
import zipfile


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Rename SMPL-X files to be loadable by smplx')
    
    parser.add_argument('smplx_zip', type=str, default='/path/to/smplx_lockedhead_20230207.zip')
    
    args = parser.parse_args()
    
    model_folder = 'models'
    smplx_zip = args.smplx_zip
    
    # Unzip the SMPL model
    print('Unzipping SMPL model...')
    with zipfile.ZipFile(smplx_zip, 'r') as zip_ref:
        zip_ref.extractall(model_folder)

    os.makedirs(os.path.join(model_folder, 'smplx'), exist_ok=True)
    
    # Rename the SMPL-X files
    print('Renaming SMPL-X files...')
    shutil.move(os.path.join(model_folder, 'models_lockedhead/smplx/SMPLX_FEMALE.npz'),
                os.path.join(model_folder, 'smplx', 'SMPLX_FEMALE.npz'))
    
    shutil.move(os.path.join(model_folder, 'models_lockedhead/smplx/SMPLX_MALE.npz'),
            os.path.join(model_folder, 'smplx', 'SMPLX_MALE.npz'))
    
    shutil.move(os.path.join(model_folder, 'models_lockedhead/smplx/SMPLX_NEUTRAL.npz'),
        os.path.join(model_folder, 'smplx', 'SMPLX_NEUTRAL.npz'))
    
    print("Cleaning up...")
    shutil.rmtree(os.path.join(model_folder, 'models_lockedhead'))
    print("Done!")
    print(" The SMPL-X models were placed in models/. You can now use the smplx library to load the SMPL-X model.")