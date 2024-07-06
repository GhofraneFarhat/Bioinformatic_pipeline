from setuptools import setup, find_packages


with open ("pip\README.md", "r") as f:
    long_description = f.read()

setup(
    name='plaspipe',
    version='0.1',
    packages=find_packages(where="pip"),
    package_dir={'': 'pip'},
    license="MIT",
    install_requires=[
        'networkx>=2.7',
        'pandas>=2.0.0',
        'matplotlib>=3.7.0',
        'seaborn>=0.12.2',
        'biopython>=1.81',
        'scipy>=1.10.1',
        'numpy>=1.24.2',
        'scikit-learn>=0.23.1',
        'tensorflow>=2.8.0',
        'spektral>=1.0.8',
        'gurobipy>=9.1.2'
    ],
    package_data={
        'plaspipe': ['../test/*.yaml']
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'plaspipe=src.plaspipe:main',
        ],
    },
    author='Ghofrane Farhat',
    author_email='ghofrane_farhat@sfu.ca',
    description='Bioinformatic pipeline for plasmid analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires=">=3.10",
    url='https://github.com/GhofraneFarhat/Bioinformatic_pipeline',
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
)
