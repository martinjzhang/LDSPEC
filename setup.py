import setuptools

exec(open("gdreg/version.py").read())

setuptools.setup(
    name="gdreg",
    version=__version__,
    description="Gene-level directional effect regression",
    url="https://github.com/martinjzhang/GDReg",
    author="Martin Jinye Zhang",
    author_email="martinjzhang@gmail.com",
    license="MIT",
    packages=["gdreg"],
    zip_safe=False,
)
