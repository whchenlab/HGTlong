"""Compatibility metadata for older HPC setuptools installations."""
from setuptools import setup

setup(
    name='HGTlong',
    version='0.1.3rc3',
    author='Weihua Chen',
    author_email='weihuachen@hust.edu.cn',
    license='CC-BY-NC-4.0',
    url='https://github.com/whchenlab/HGTlong',
    description='Read-resolved candidate HGT loci from ONT molecules',
    packages=['hgtlong'],
    python_requires='>=3.9',
    entry_points={'console_scripts':['hgtlong=hgtlong.cli:main']},
)
