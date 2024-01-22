# SMPL to AddBiomechanics

This repository contains Python code for generating input data for [AddBiomechanics](https://addbiomechanics.org/) from SMPL sequences. [AddBiomechanics](https://addbiomechanics.org/) is a free online tool designed to align an OpenSim Skeleton model to motion capture data.

# Installation

First install the dependancies and the package with the following line:
```pip install -e .```


## SMPL

Download the file: SMPL_python_v.1.1.0.zip from the [SMPL download page](https://smpl.is.tue.mpg.de/). And run:

```
cd ../SMPL2AddBiomechanics
python scripts/setup_smpl.py /path/to/SMPL_python_v.1.1.0.zip  
```

## BSM 

Download the file: skel_models_v1.x.zip from the [SKEL download page](https://skel.is.tue.mpg.de/). And run:

```
cd ../SMPL2AddBiomechanics
python scripts/setup_bsm.py /path/to/skel_models_v1.1.zip
```


## OSSO (Optional)
To generate a personalized custom markers on the OpenSim model, you will need [OSSO](https://github.com/MarilynKeller/OSSO/). Please check the [OSSO installation instructions](https://github.com/MarilynKeller/OSSO/blob/main/installation.md) to set it up.


# Usage

Below is an example that takes the SMPL motion `models/bsm/sample_motion/01_01_poses.npz` and generates measurements and synthetic motion capture data from the sequence.

```bash
python smpl2ab/smpl2addbio.py -i  models/bsm/sample_motion/01 
``` 

This generates an output folder containing the following files:
```
output
└── 01
    ├── _subject.json
    └── trials
        └── 01_01_poses.trc
```

The `_subject.json` file contains the estimated subject measurements :
`{"sex": "male", "massKg": "68.94", "heightM": "1.80"}`
, and the `trc` file contains the synthetic motion capture data for the input sequences.

Those informations can be used as input for [AddBiomechanics](https://dev-addbiomechanics.org/) to align the OpenSim model to the motion capture data.
Note that you will need to use the developper version of https://dev-addbiomechanics.org/ that support the BSM model.

In the 'Upload Custom OpenSim Model' section, upload the BSM model `models/bsm/bsm.osim`.

In "Motion Capture Files" section, upload the `trc` file from `output/01/trials/01_01_poses.trc`.


## Custom markers

The markers location wrt the bones depend on the body shape. While AddBiomechanics does optimize the markers location, we propose a method that leverages OSSO to generate a personalized custom markers set on the OpenSim model. 
To use this option, add the '--osso' flag to the command line. This will generate a `bsm.osim` file in the output folder that contains a version of the BSM model with personalized markers. 

```bash 
python smpl2ab/smpl2addbio.py -i  models/bsm/sample_motion/01 --osso --display
```

Then proceed as above to upload the data to [AddBiomechanics](https://dev-addbiomechanics.org/). But this time, upload the `bsm.osim` file from `output/01/bsm.osim` in the "Upload Custom OpenSim Model" section.

## Running AddBiomechanics locally

If you want to run AddBiomechanics locally, you can download and install the source code from the [AddBiomechanics github page](https://github.com/keenon/AddBiomechanics).

The script [AddBiomechanics/blob/main/server/engine/src/engine.py](https://github.com/keenon/AddBiomechanics/blob/main/server/engine/src/engine.py) runs the optimization given a folder of data.

# Running the code on your own data

You will need
- A SMPL motion sequence
- A biomechanical skeleton model .osim with markers defined on it (BSM will be used by defaut)
- A dictionary giving, for each marker of the .osim, the corresponding vertex index on SMPL (The BSM correspondance will be used by default)

## SMPL sequences
You can download SMPL sequences from the [AMASS](https://amass.is.tue.mpg.de/) download page, by clicking the `SMPL+H` button. 
Here we show an example for the subject `10` of `CMU`. The AMASS data consists in a folder with this structure:
```
CMU
└── 10
    ├── 10_01_poses.npz
    ├── 10_02_poses.npz
    ├── 10_03_poses.npz
    ├── 10_04_poses.npz
    └── 10_05_poses.npz
```

Each .npz is a sequence.

You should pass the path to the subject folder to the `smpl2addbio.py` script.


## Biomechanical model

By default, you can use our BSM model:
Download bsm.osim and the Geometry folder from the [SKEL project page](https://skel.is.tue.mpg.de/) `download page - > "Download Models" button`, and place it in the 'models' folder.

You can also use your own OpenSim model:

python smpl2ab/smpl2addbio.py -i  models/bsm/sample_motion/01 --osso --display

## Markers on SMPL

For each marker of the .osim model, you need to provide the corresponding vertex index on SMPL.
By default, we use the BSM marker set. You can find it as a dictionnary in `SMPL2AddBiomechanics/smpl2ab/data/bsm_markers.yaml`.

If your model has a different marker set, you can create your own dictionnary similar `SMPL2AddBiomechanics/smpl2ab/data/bsm_markers.yaml`.
This dictionnary must contain the same markers as your .osim model. To help you to create this dictionnary, you can run the following script:
```
python smpl2ab/utils/kin_helpers.py -o  /path/to/your/model.osim  -m -mr -D
```

This will print all the markers of your OpenSim model,the bone they are rigged to and display each markers on the skeleton mesh with their label. Note that to display the model, you will need to install [aitviewer-skel](https://github.com/MarilynKeller/aitviewer-skel).


# Citation
If you use this software, please cite the following work and software:

```
@inproceedings{keller2023skel,
  title = {From Skin to Skeleton: Towards Biomechanically Accurate 3D Digital Humans},
  author = {Keller, Marilyn and Werling, Keenon and Shin, Soyong and Delp, Scott and 
            Pujades, Sergi and C. Karen, Liu and Black, Michael J.},
  booktitle = {ACM ToG, Proc.~SIGGRAPH Asia},
  volume = {42},
  number = {6},
  month = dec,
  year = {2023},
}
```

## License

This code and model are available for non-commercial scientific research purposes as defined in the [LICENSE.txt](LICENSE.txt) file.


# Contact 

For any question about SKEL loading, please contact skel@tuebingen.mpg.de.

For commercial licensing, please contact ps-licensing@tue.mpg.de