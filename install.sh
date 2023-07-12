#!/bin/bash

## build environment

virtualenv --python python3 venv

source venv/bin/activate

pip install --upgrade pip

pip install -r packages.txt

deactivate
