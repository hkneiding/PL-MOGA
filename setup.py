from setuptools import setup, find_packages

setup(
    name='gapy',
    version='1.0',

    url='https://github.com/hkneiding/PL-MOGA',
    author='Hannes Kneiding',
    author_email='hannes.kneiding@outlook.com',
    packages=find_packages(exclude=['test'])
)
