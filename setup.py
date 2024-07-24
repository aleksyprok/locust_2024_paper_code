"""
This is the setup script for the iaea_fec_code package.

To install this package using pip, run the following command:
    pip install -e .

After installation, other modules can import this package for testing purposes.
"""

from setuptools import setup, find_packages

setup(
    name='nuclear_fusion_code',
    version='0.1.0',
    packages=find_packages(),
)
