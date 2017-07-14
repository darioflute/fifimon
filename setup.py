#!/usr/bin/env python

from distutils.core import setup

setup(name='fifimon',
      version='0.6.alpha',
      description='FIFI-LS Monitor',
      long_description='The package monitors FIFI-LS data',
      author='Dario Fadda',
      author_email='darioflute@gmail.com',
      url='https://github.com/darioflute/fifimon.git',
      download_url='https://github.com/darioflute/fifimon',
      license='GPLv3+',
      packages=['fifimon'],
      scripts=['bin/fifimon']
     )
