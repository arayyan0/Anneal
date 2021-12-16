./clean_jobs.sh
vim sweeps/jtau_sweep.py
python3 sweeps/jtau_sweep.py
tail -n1 jobrun_*.sh | bash
clear
