from setuptools import setup, find_packages

setup(
    name="CASSIA",
    version="0.1.8",
    packages=find_packages(),
    package_data={
        'CASSIA': ['data/*.csv'],  # Include all CSV files in data directory
    },
    include_package_data=True,  # This tells setuptools to include package_data
    install_requires=[
        "numpy",
        "pandas",
        "openai>=1.0.0",
        "anthropic>=0.3.0",
        "requests>=2.25.0",
    ],
    author="Elliot Yixuan Xie",
    author_email="xie227@wisc.edu",
    description="CASSIA (Cell type Annotation using Specialized System with Integrated AI) is a Python package for automated cell type annotation in single-cell RNA sequencing data using large language models.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/elliotxe/CASSIA",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
)

