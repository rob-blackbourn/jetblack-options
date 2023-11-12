from setuptools import setup

setup(
    name="jetblack-options",
    description="Reference implementations for option pricing formula",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    version="0.1.0",
    url="https://github.com/rob-blackbourn/jetblack-options",
    author="Rob Blackbourn",
    author_email="rob.blackbourn@gmail.com",
    license="Apache-2.0",
    classifiers=[
        "Operating System :: OS Independent",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.8",
    ],
    packages=["jetblack_options"],
    install_requires=[],
)
