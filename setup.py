import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="starcatpy",
    version="1.0.6",
    author="Dylan Kotliar, Michelle Curtis",
    author_email="dylkot@gmail.com, curtism@broadinstitute.org",
    description="Implements *CellAnnotator (aka *CAT/starCAT), annotating scRNA-Seq with predefined gene expression programs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/immunogenomics/starCAT",
    project_urls={
        "Bug Tracker": "https://github.com/immunogenomics/starCAT/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    include_package_data=True,
    package_data={
        "starcat": ["current_references.tsv"],
    },
    entry_points={
        'console_scripts': [
            'starcat = starcat:main',
        ],
    },
    install_requires=[
   'scikit-learn>=1.0',
   'anndata',
   'pandas',
   'numpy',
   'scipy',
   'pyyaml',
   'requests'
   ]
)
