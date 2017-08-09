#!/usr/bin/env python

from distutils.core import setup
import setuptools

setup(name='fifimon',
      version='0.22.alpha',
      description='FIFI-LS Monitor',
      long_description='The package monitors FIFI-LS data',
      author='Dario Fadda',
      author_email='darioflute@gmail.com',
      url='https://github.com/darioflute/fifimon.git',
      download_url='https://github.com/darioflute/fifimon',
      license='GPLv3+',
      packages=['fifimon'],
      scripts=['bin/fifimon'],
      include_package_data=True,
      package_data={'fifimon':['icons/*.jpg','copyright.txt','CalibrationResults.csv']}
     )
