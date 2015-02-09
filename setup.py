from setuptools import setup, find_packages

setup(
    name = "DumbView",
    version="0.3.1",
    author="Pacific Biosciences",
    author_email="dalexander@pacificbiosciences.com",
    scripts = ["bin/dumbview", "bin/dumbview-ccs"],
    packages = find_packages(),
    zip_safe = False,
    install_requires=[
        "GenomicConsensus >= 0.9.1",
        "pbcore >= 0.9.3",
        "numpy >= 1.6.0",
        "h5py >= 2.0.1"
        ]
    )
