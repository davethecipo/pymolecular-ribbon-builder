#!/usr/bin/env python

from setuptools import setup

setup(name='molecular-ribbon-builder',
      version='0.1',
      description='build molecular ribbons starting from oligomers',
      author='Davide Olianas',
      author_email='ubuntupk@gmail.com',
      url='',
      download_url = '',
      keywords = ['molecule', 'ribbon', '2D'],
      packages=['ribbonbuilder'],
      entry_points={
          'console_scripts': [
              'molribbon2d = ribbonbuilder.molribbon2d:main']
      },
      install_requires=['PyYAML', 'numpy', 'openbabel'],
      extras_require={
          'dev': [ 'Sphinx', 'alabaster', 'pytest' ]
      }
     )
