from setuptools import setup, find_packages

REQUIREMENTS = ['scipy', 'numpy','pandas', 'sklearn']

setup(
    name='oloPATH',
    version='1.0',
    description='Package for metabolomic pathways analysis',
    author='Neus Pou Amengual',
    author_email='np@olobion.ai',
    url='https://github.com/oloBion/oloPATH',
    install_requires = REQUIREMENTS,
    packages = find_packages(),
    package_data={
    'olopath': [
        'data/ChEBI/*.json.zip',
        'data/reactome/*.json.zip',
        'data/metabolic_pathways/*.json.zip',
        ],
    },
    python_requires='>=3.10'
    )