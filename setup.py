"""
GeoVac: Geometric Vacuum Quantum Solver
A sparse-matrix quantum chemistry solver using AdS5 paraboloid lattice discretization.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name='geovac',
    version='0.7.0',
    author='J. Loutey',
    author_email='jloutey@example.com',  # Update with your actual email
    description='Topological Hartree-Fock quantum solver with universal constant -1/16',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jloutey-hash/geovac',
    license='MIT',
    
    packages=find_packages(exclude=['tests', 'old_research_archive', 'geometric-solver']),
    
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'networkx>=2.6.0',
    ],
    
    extras_require={
        'dev': [
            'pytest>=6.0',
            'matplotlib>=3.3.0',
        ],
        'visualization': [
            'matplotlib>=3.3.0',
        ],
    },
    
    python_requires='>=3.8',
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
    ],
    
    keywords=[
        'quantum-computing',
        'holographic-principle',
        'sparse-matrix',
        'physics',
        'quantum-chemistry',
        'computational-physics',
        'ads-cft',
        'lattice-theory',
    ],
    
    project_urls={
        'Documentation': 'https://github.com/jloutey-hash/geovac/wiki',
        'Source': 'https://github.com/jloutey-hash/geovac',
        'Bug Reports': 'https://github.com/jloutey-hash/geovac/issues',
    },
)
