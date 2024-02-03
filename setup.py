import setuptools

setuptools.setup(
    name="toolkit",
    version="1.0",
    packages=["toolkit", "weights"],
    install_requires=["numpy", "scipy", "healpy", "pixell", "matplotlib", "astropy"],
)
