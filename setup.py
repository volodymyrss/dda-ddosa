#!/usr/bin/env python

from setuptools import setup

setup(name='dda-ddosa',
      version='1.0.1',
      description='ddosa',
      author='V.S.',
      py_modules=['ddosa'],
      install_requires=[
         "pilton>=1.0.1",
         "integral-data-mirror>=1.1.0",
      ]
     )
