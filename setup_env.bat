@echo off
echo Setting up Python environment for DNA Sequencing project...

:: Create virtual environment
python -m venv venv

:: Activate virtual environment
call venv\Scripts\activate

:: Upgrade pip
python -m pip install --upgrade pip

:: Install requirements
pip install -r requirements.txt

:: Install the project in development mode
pip install -e .

:: Create project structure
mkdir src\generators src\algorithms src\utils src\config
mkdir tests\test_generators tests\test_algorithms tests\test_utils
mkdir data\raw data\results
mkdir docs\report
mkdir notebooks

:: Create __init__.py files
type nul > src\__init__.py
type nul > src\generators\__init__.py
type nul > src\algorithms\__init__.py
type nul > src\utils\__init__.py
type nul > src\config\__init__.py
type nul > tests\__init__.py

echo Environment setup complete!
echo To activate the environment, run: venv\Scripts\activate 