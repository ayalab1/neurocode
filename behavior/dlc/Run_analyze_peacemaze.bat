@echo on
call C:\ProgramData\Anaconda3\Scripts\activate.bat
call conda activate C:\Users\Cornell\.conda\envs\DEEPLABCUT
python "D:\github\neurocode\behavior\dlc\analyze_peacemaze.py"

