"""
++ mkalgo
++ saif.778@gmail.com
"""

import re
import ast
from setuptools import setup

install_requires = [
    'numpy',
]

tests_require = []

setup(
    name='mkalgo',
    version='0.2',
    url='https://github.com/saifuddin778/mkalgo/',
    include_package_data=True,
    description='mkalgo',
    author='saifuddin778',
    install_requires=install_requires,
    extras_require={'test': tests_require},
    packages=['mkalgo'],
    entry_points='''''',
    tests_require=tests_require,
)
