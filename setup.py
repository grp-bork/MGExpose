# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys

from mgexpose import __version__ as mgexpose_version
here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="mgexpose"
	version = mgexpose_version

	if sys.version_info.major != 3:
		raise EnvironmentError(f"""{name} requires Python3, and is not compatible with Python2.""")

	setup(
		name=name,
		version=version,
		description=description,
		long_description=long_description,
		url="https://github.com/cschu/mgexpose",
		author="Christian Schudoma, Anastasia Grekova, Supriya Khedkar",
		author_email="christian.schudoma@embl.de, anastasiia.grekova@embl.de, khedkar@bioquant.uni-heidelberg.de ",
		license="MIT",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific/Engineering :: Bio-Informatics",
			"License :: OSI Approved :: MIT License",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Programming Language :: Python :: 3.10",
			"Programming Language :: Python :: 3.11",
			"Programming Language :: Python :: 3.12",
			"Programming Language :: Python :: 3.13",
		],
		zip_safe=False,
		keywords="microbial genome mobile genetic element detection",
		packages=find_packages(exclude=["test"]),
		entry_points={
			"console_scripts": [
				"mgexpose=mgexpose.__main__:main",				
			],
		},
		package_data={},
		include_package_data=True,
		data_files=[],
		install_requires=[
			"pyhmmer",
			"pyrodigal",
		],
	)
