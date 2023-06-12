import setuptools

exec(open("ldspec/version.py").read())

setuptools.setup(
    name="ldspec",
    version=__version__,
    description="LD-score SNP-pair effect correlation regression",
    url="https://github.com/martinjzhang/LDSPEC",
    author="Martin Jinye Zhang",
    author_email="martinjzhang@gmail.com",
    license="MIT",
    packages=["ldspec"],
    zip_safe=False,
)
