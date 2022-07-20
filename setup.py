#!/usr/bin/env python

from setuptools import find_packages, setup

classifiers = [
    'Development Status :: 1 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'VXdetector'


setup(name='vxdetector',
      version='1.0.0',
      license='BSD',
      #description=description,
      #long_description=long_description,
      #keywords=keywords,
      classifiers=classifiers,
      author="Johannes Groos",
      author_email="Johannes.Groos@bio.uni-giessen.de",
      maintainer="http://jlab.bio",
      maintainer_email="stefan.m.janssen@gmail.com",
      url='https://github.com/jlab/algorithm_vxdetector/',
      test_suite='nose.collector',
      packages=find_packages(),
      install_requires=[
          #'click >= 6',
          #'scikit-bio >= 0.4.0',
      ],
      #extras_require={'test': ["nose", "pep8", "flake8"],
      #                      'coverage': ["coverage"]}
)
