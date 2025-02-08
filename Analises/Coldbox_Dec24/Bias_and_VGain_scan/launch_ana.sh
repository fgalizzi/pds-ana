MODULE=$1
echo -e "\n \n Analysing module: $MODULE  \n \n"

eval 'root -e ".L ./VGainScans_ana.cpp+" -e "auto a = cla()" -e "VGainScans_ana(a, \"config/config_module_$MODULE.json\")" '
