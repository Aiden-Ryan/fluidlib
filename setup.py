from setuptools import setup, find_packages

setup(
    name = 'fluidlib',
    version = '0.1.0',
    description='fluid library with transient thermal solver',
    author = 'Aiden Ryan',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy', 
        'matplotlib',
    ]  
)