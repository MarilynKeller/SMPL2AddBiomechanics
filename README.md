# SMPL to AddBiomechanics

This repository contains Python code for generating input data for [AddBiomechanics](https://addbiomechanics.org/) from SMPL sequences. [AddBiomechanics](https://addbiomechanics.org/) is a free online tool designed to align an OpenSim Skeleton model to motion capture data.

# Installation

First install the dependancies and the package with the following line:
```pip install -e .```

# OSSO
To generate a personalized custom markers on the OpenSim model, you will need [OSSO](https://github.com/MarilynKeller/OSSO/). Please check the [OSSO installation instructions](https://github.com/MarilynKeller/OSSO/blob/main/installation.md) to set it up.

## Usage

## Example

```bash
python main.py -i /path/to/smpl/sequences -o /path/to/output/folder 
```

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