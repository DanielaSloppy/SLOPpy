from setuptools import setup

# Inspired from here:
# https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
# https://realpython.com/pypi-publish-python-package/#prepare-your-package-for-publication

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="SLOPpy-package",
    version='1.5.1',
    author="Daniela Sicilia, Luca Malavolta, et al.",
	author_email = 'daniela.sicilia@inaf.it, luca.malavolta@unipd.it',
	url = 'https://github.com/LucaMalavolta/SLOPpy',
	packages =['SLOPpy', 'SLOPpy.subroutines', 'SLOPpy.instruments'],
	license = 'MIT License',
	description ='SLOPpy: Spectral Lines Of Planets with python',
    long_description=long_description,
    long_description_content_type='text/markdown',
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
		'Programming Language :: Python :: 3'
		],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'SLOPpy_run=SLOPpy.sloppy_run:sloppy_run',
        ]
    },
    zip_safe=False,
    install_requires=[
        'numpy==1.26.4',
        'numba>=0.55.1',
        'scipy>=1.8.1',
        'matplotlib>=3.5.2',
        'astropy>=5.1',
        'astroquery>=0.4',
        'pyerfa>=2.0',
        'argparse>=1.4',
        'oyaml>=1.0',
        'emcee>=3.1.2',
        'pyyaml',
        'h5py>=3.7.0',
        'tqdm>=4.60',
        'pygtc>=0.4.1',
        'jaxlib==0.4.31',
        'jax==0.4.31',
        'tinygp>=0.3.0',
        'PyAstronomy>=0.18',
        'ordered-set',
        'scikit-learn'
    ],
    setup_requires=['setuptools']
)
