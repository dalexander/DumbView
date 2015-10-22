from setuptools import setup, find_packages

globals = {}
execfile("DumbView/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]

setup(
    name = "DumbView",
    version=__VERSION__,
    author="Pacific Biosciences",
    author_email="dalexander@pacificbiosciences.com",
    scripts = ["bin/dumbview", "bin/dumbview-ccs"],
    packages = find_packages(),
    zip_safe = False,
    install_requires=[
        "GenomicConsensus >= 0.9.1",
        "pbcore >= 1.2.3",
        "numpy >= 1.6.0",
        "h5py >= 2.0.1"
        ]
    )
