MODULE=$1
echo -e "\n \n Analysing module: $MODULE  \n \n"

eval 'root -e ".L ./Bin.cpp+" -e "auto a = cla()" -e "Bin(a, \"config/config_module_$MODULE.json\")" '

