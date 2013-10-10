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
        "pbcore >= 0.8.0",
        "numpy >= 1.6.0",
        "h5py >= 2.0.1"
        ]
    )
