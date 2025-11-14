# amp_desige

### How to run it on your own machine

1. # Install the python version don't newer then v3.11

    https://www.python.org/downloads/

2. # active virtual environment

   ```
   $ python -m venv amp_env
   $ Set-ExecutionPolicy -Scope Currentuser -ExecutionPolicy Remotesigned
   $ .\amp_env\Scripts\Activate.ps1
   ```

3. # Install the requirements

   ```
   $ pip install wheel setuptools
   $ pip install -r requirements.txt
   ```

4. # Run the app

   ```
   $ streamlit run amp_design_app.py
   ```
