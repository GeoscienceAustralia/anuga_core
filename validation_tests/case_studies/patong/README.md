### TO RUN THIS CASE STUDY (16/01/2015):

```
python data_download.py
mpirun -np XX python run_model.py > outfile.log & ## Here XX = number of processors
```

This runs a Patong simulation to produce an SWW file. 

### 6/06/2013 GD modifications

- Updated names of functions to reflect code refactoring

- Run 'python test_patong_scenario.py' to get data. 
    I could successfully get the data.
    The script then tried to the run the models, but it failed.

- I then ran 
```
python build_elevation.py
python build_urs_boundary.py
```

- Set all 'use_cache' options to 'False', as it was causing errors.

- `test_patong_scenario.py` was still having problems running the models. However,
    I found I could run `python run_model.py` successfully.
    Thus, moved `test_patong_scenarios.py` to `download_data_for_patong_scenarios.py`,
    and cut out the code that tries to do the testing. 
    For the purposes of the validation document, this is probably a reasonable idea.

- Debugged it, and got rid of `setup_model.py` to simplify the code

- Got it running in parallel

### 16/01/2015 SR modifications

- Now the data lives in an svn repository which you download via the script data_download.py
