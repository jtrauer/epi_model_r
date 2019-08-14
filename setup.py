from setuptools import setup, find_packages

setup(
    name='summer',
    version='1.0.3',
    packages=["summer"],
    url='https://github.com/jtrauer/summer',
    license='MIT',
    author='James Trauer',
    author_email='james.trauer@monash.edu',
    install_requires=['scipy==1.1.0',
                      'graphviz>=0.4.10',
                      'SQLAlchemy>=1.1.18',
	              #'pymc==2.3.6',	
                      'pandas==0.24.2',],
    description='General structure for creating epidemiological models in R'
)
