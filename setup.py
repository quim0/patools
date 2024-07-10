import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="patools",
    version="0.1.0",
    author="Quim Aguado",
    author_email="",
    description="Tools for working with pairwise sequence alignment files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/quim0/patools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    setup_requires=['wheel'],
    install_requires=[
        'numpy',
        'matplotlib',
        'rich',
        'tqdm'
    ],
    entry_points = {
        'console_scripts': [
            'catcigar=catcigar.catcigar:catcigar',
            'checkalign=checkalign.checkalign:checkalign'],
    }
)