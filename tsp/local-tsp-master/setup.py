# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='tsp_local',
    version='0.1.0',
    description='Implementation of TSP heuristics',
    long_description=readme,
    author='Arthur Mahéo',
    author_email='arthur.maheo@anu.edu.au',
    url='https://gitlab.com/Soha/local-tsp',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')))
