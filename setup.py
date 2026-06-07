from setuptools import setup, find_packages

setup(
    name="UfukMSA-221201022",
    version="0.1.0",
    author="Ufuk Dogan",
    author_email="221201022@rumeli.edu.tr",
    description="MAFFT-based Multiple Sequence Alignment",
    packages=find_packages(),
    install_requires=["numpy"],
    python_requires=">=3.7",
)
