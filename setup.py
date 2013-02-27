from setuptools import setup, find_packages

setup(
    name = "DumbView",
    version="0.1.0",
    author="Pacific Biosciences",
    author_email="dalexander@pacificbiosciences.com",
    scripts = ["bin/dumbview"],
    packages = find_packages(),
    zip_safe = False,
    install_requires=[
        "pbcore >= 0.5.0",
        "numpy >= 1.6.0",
        "h5py >= 1.3.0"
        ]
    )
