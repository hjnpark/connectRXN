from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='connectRXN',
    description='Connects unit reactions to create unique pathways',
    url='https://github.com/hjnpark/connectRXN',
    author='Heejune Park',
    packages=find_packages(),
    package_data={'': ['*.ini']},
    include_package_data=True,
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={'console_scripts': [
        'connect-rxn = connect.connect:main',
    ]},
    install_requires=[
        'numpy>=1.20',
        'networkx',
        'matplotlib',
        'pyvis',
        'rdkit',
        'pillow',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
    zip_safe=True,
)
