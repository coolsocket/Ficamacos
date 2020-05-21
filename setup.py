import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    'numpy>=1.11.1',
    'fft'
]

setuptools.setup(
    name="Ficamacos", # Replace with your own username
    version="0.1.0",
    author="Yehao Liu",
    author_email="594310297@qq.com",
    description="Option pricing with COS, Carr Madan and Filon's method",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/coolsocket/Ficamacos",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7'
)
