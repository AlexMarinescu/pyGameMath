import os
from setuptools import setup, find_packages
setup(
  name = 'pgm',
  packages = find_packages(where="pgm"),
  package_dir = {'':'pgm'},
  version = '0.1.0',
  description = 'Math library for game programming in python. ',
  author = 'Alex Marinescu',
  author_email = 'ale632007@gmail.com',
  license='BSD 2-Clause',
  url = 'https://github.com/explosiveduck/pyGameMath',
  download_url = 'https://github.com/explosiveduck/pyGameMath/tarball/v0.1.0',
  install_requires=['six'],
  keywords = ['math', 'game', 'library'],
  classifiers = [
            'Development Status :: 0.1.0',
            'Intended Audience :: Developers',
            'Topic :: Libraries :: Graphics :: Mathematics :: Software Development',
            'License :: OSI Approved :: BSD 2-Clause "Simplified" or "FreeBSD" license (BSD-2-Clause)',

            'Programming Language :: Python',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            ],
)