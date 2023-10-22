from setuptools import setup

setup(
    name='ww3tools',
    version='1.2.0',
    description='Python tools and utilities for WAVEWATCHIII post-processing and validation.',
    authors='Ricardo M. Campos, Ali Abdolali, Matthew Masarik, Saeideh Banihashemi, Maryam Mohammadpour, Jessica Meixner', 
    author_email='ricardo.campos@noaa.gov',
    author_website='https://www.aoml.noaa.gov/people/ricardo-campos/',
    credits='National Oceanic and Atmospheric Administration (NOAA) and Cooperative Institute For Marine And Atmospheric Studies (CIMAS)',
    url='https://github.com/NOAA-EMC/WW3-tools',
    license='LGPLv3',
    packages=["ww3tools",
        "ww3tools.downloadobs"],
    install_requires=[
        "numpy",
        "scipy",
        "statistics",
        "matplotlib",
        "cartopy",
        "basemap",
        "netCDF4",
        "pandas",
        "xarray",
        "pyresample",
        "salem",
        "regionmask",
        "geopy",
        "pathlib2",
        "tqdm",
        "yaml",
        "boto3"]
)

