from setuptools import setup, find_packages

setup(
    name = "DumbView",
    version="0.2.0",
    author="Pacific Biosciences",
    author_email="dalexander@pacificbiosciences.com",
    scripts = ["bin/dumbview"],
    packages = find_packages(),
    zip_safe = False,
    install_requires=[
        "ConsensusCore >= 0.7.6",
        "GenomicConsensus >= 0.8.0",
        "pbcore >= 0.8.0",
        "numpy >= 1.6.0",
        "h5py >= 2.0.1"
        ]
    )
