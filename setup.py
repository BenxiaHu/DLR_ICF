from setuptools import setup
from setuptools import find_packages

version_py = "DLR_ICF/_version.py"
exec(open(version_py).read())

setup(
    name="DLR_ICF", # Replace with your own username
    version=__version__,
    author="Benxia Hu",
    author_email="hubenxia@gmail.com",
    description="DLR and ICF analysis",
    long_description="analyze DLR and ICF",
    url="https://pypi.org/project/DLR_ICF/",
    entry_points = {
        "console_scripts": ['DLR_ICF_main = DLR_ICF.DLR_ICF_main:main',
                            'DLR_ICF_comparison = DLR_ICF.DLR_ICF_comparison:main',
                            'DLR_ICF_separation = DLR_ICF.DLR_ICF_separation:main',
                            'ICF_chromatin = DLR_ICF.ICF_chromatin:main',]
        },
    python_requires = '>=3.6',
    packages = ['DLR_ICF'],
    install_requires = [
        'numpy',
        'pandas',
        'argparse',
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    zip_safe = False,
  )