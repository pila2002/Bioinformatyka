from setuptools import setup, find_packages

setup(
    name="dna-sequencing",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        # Will be read from requirements.txt
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="DNA Sequencing by Hybridization with Metaheuristics",
    python_requires=">=3.8",
) 