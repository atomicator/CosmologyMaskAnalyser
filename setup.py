import setuptools

setuptools.setup(
    name="toolkit",
    version="1.0",
    description="This should make conda talk to the file",
    author="me",
    packages=["toolkit"],
    install_requires=["numpy", "scipy", "healpy", "pixell", "matplotlib", "astropy"],
)
