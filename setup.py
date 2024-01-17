from setuptools import setup, find_packages

requirements =[
    "numpy",
    "torch",
    "smplx",
    "nimblephysics==0.10.32",
    "tqdm",
    "rtree",
]

setup(
    name="smpl2ab",
    version="0.1",
    description=("Generating input data for [AddBiomechanics](https://addbiomechanics.org/) from SMPL sequences."),
    author="Marilyn Keller",
    author_email="marilyn.keller@tuebingen.mpg.de",
    packages=find_packages(),
    install_requires=requirements,
)

