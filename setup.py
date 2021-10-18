"""Titan aerosols setup."""

from pathlib import Path

from setuptools import setup


HERE = Path(__file__).parent
README = (HERE / 'README.md').read_text()


setup(
    name='titan-aerosols',
    version='0.3.0',
    description='Titan aerosols models',
    author='Benoit Seignovert, Pascal Rannou, Loic Rossi',
    author_email='python@seignovert.fr',
    url='http://github.com/seignovert/python-titan-aerosols',
    license='MIT',
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20',
    ],
    packages=['aerosols'],
    include_package_data=True,
    keywords=['Titan', 'Aerosols', 'Scattering models', 'Mie', 'Fractal'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],
    long_description=README,
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'fractal_tholins=aerosols.cli:cli_fractal_tholins'
        ]
    },
)
