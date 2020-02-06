import os
from setuptools import setup, find_packages

# Read package installation requirements from 'requirements.txt'
file_dir = os.path.dirname(os.path.abspath(__file__))
reqs_path = os.path.join(file_dir, "requirements.txt")
with open(reqs_path) as f:
    reqs_text = f.read()
    package_reqs = [
        line for line in reqs_text.split("\n") if line and not line.strip().startswith("#")
    ]

setup(
    name="summer_py",
    version="1.1.0",
    packages=find_packages("."),
    url="https://github.com/jtrauer/summer_py",
    license="MIT",
    author="James Trauer",
    author_email="james.trauer@monash.edu",
    install_requires=package_reqs,
    description="General structure for creating epidemiological models in Python",
)
