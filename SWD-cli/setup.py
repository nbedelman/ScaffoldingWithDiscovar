"""Packaging settings."""


from codecs import open
from os.path import abspath, dirname, join
from subprocess import call

from setuptools import Command, find_packages, setup

from SWD import __version__


this_dir = abspath(dirname(__file__))
with open(join(this_dir, 'README.rst'), encoding='utf-8') as file:
    long_description = file.read()


class RunTests(Command):
    """Run all tests."""
    description = 'run tests'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run all tests!"""
        errno = call(['py.test', '--cov=SWD', '--cov-report=term-missing'])
        raise SystemExit(errno)


setup(
    name = 'SWD',
    version = __version__,
    description = 'A command line program in Python intended to use DISCOVAR-type data to order and orient reference scaffolds.',
    long_description = 'A command line program in Python intended to use DISCOVAR-type data to order and orient reference scaffolds.',
    url = 'https://github.com/nbedelman/ScaffoldingWithDiscovar',
    author = 'Nate Edelman',
    author_email = 'nedelman@g.harvard.edu',
    license = 'UNLICENSE',
    classifiers = [
        'Intended Audience :: Scientists',
        'Topic :: Utilities',
        'License :: Public Domain',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords = 'cli',
    packages = find_packages(exclude=['docs', 'tests*']),
    install_requires = ['docopt'],
    extras_require = {
        'test': ['coverage', 'pytest', 'pytest-cov'],
    },
    entry_points = {
        'console_scripts': [
            'SWD=SWD.cli:main',
        ],
    },
    cmdclass = {'test': RunTests},
)
