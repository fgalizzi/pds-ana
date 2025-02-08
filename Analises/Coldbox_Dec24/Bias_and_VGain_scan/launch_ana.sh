MODULE=$1
echo -e "\n \n Analysing module: $MODULE  \n \n"

eval 'root  -e ".L ./VGainScans_ana.cpp+" -e "#include \"config/class_configs_module_$MODULE.hpp\" " -e "#include \"config/module_$MODULE.hpp\" " -e "auto a = cla()"  '
# eval 'root -e "#include \"config/module_$MODULE.hpp\" " -e ".L ./VGainScans_ana.cpp+" -e "#include \"config/class_configs_module_$MODULE.hpp\" " -e "auto a = cla()" -e "VGainScans_ana(a)" '
