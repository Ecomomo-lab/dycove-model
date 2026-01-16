@echo off
REM Batch script to run ANUGA pre-processing and parallel simulation locally.
REM Adapted from a Slurm script for use in Windows Anaconda Prompt.

REM --- STEPS FOR RUNNING ANUGA IN PARALLEL ---
REM 1. Open Anaconda Prompt or Miniconda Prompt.
REM 2. Activate the conda environment that has Python, ANUGA, MPI (e.g., MS-MPI),
REM    mpi4py, and all other dependencies installed:
REM    conda activate your_anuga_env_name
REM 3. Navigate to the directory containing this batch file and your Python scripts:
REM    cd C:\path\to\your\anuga\project
REM 4. *** IMPORTANT: EDIT THE 'NUM_PROCESSES' VALUE BELOW! ***
REM 5. Run the script by typing its name:
REM    .\tide_channel_ANUGA_pll.bat

REM --- CONFIGURATION ---
REM Set the number of processes for MPI run below (NUM_PROCESSES).
REM Determine based on computer resources (number of cores and RAM).
REM Ensure that the ratio n_grid_cells/NUM_PROCESSES > ~1500.
REM Start conservatively (e.g., 2 or 4) and monitor TASK MANAGER (Memory usage).

SET NUM_PROCESSES=4

echo ============================================================
echo Starting Local ANUGA Workflow
echo Current time: %TIME% on %DATE%
echo Number of processes set to: %NUM_PROCESSES%
echo ============================================================
echo.

mpiexec -np %NUM_PROCESSES% python tide_channel_ANUGA_pll.py
IF %ERRORLEVEL% NEQ 0 (
    echo WARNING: Model run returned an error during execution.
    echo Check the output above and any generated log/error files.
    REM goto :eof
) else (
    echo Parallel run finished successfully according to exit code.
)
echo.

REM --- COMPLETION ---
echo ============================================================
echo ANUGA Workflow Script Finished.
echo Final time: %TIME% on %DATE%
echo Review output files and any messages above.
echo ============================================================
echo.

:eof
REM End of script
