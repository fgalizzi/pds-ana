RUN=$1
echo -e "\n \n Analysing module: $RUN  \n \n"

eval 'root -e ".L ./Calibration_ana.cpp" -e "auto a = cla()" -e "Calibration_ana(a, \"config/config_run_$RUN.json\")" '
