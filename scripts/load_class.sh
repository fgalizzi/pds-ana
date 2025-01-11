BASE_PATH="YOUR_PATH_TO/pds-ana"
echo -e "\n \n Running the general class Analyzer \n \n"

if [[ "$1" == "const" ]]; then
  eval 'root  -e "#include \"const.hpp\" " -e "std::string classe_path=\"$BASE_PATH\" " -e "#include \"$BASE_PATH/Class/_c/class_include_const.hpp\" " -e "auto a = cla()" -e "a.class_path=classe_path"'
elif [[ "$1" == "+" ]]; then
  eval 'root  -e "std::string classe_path=\"$BASE_PATH\" " -e ".L $BASE_PATH/Class/_c/class_include.hpp+" -e "#include \"test.hpp\" " -e "auto a = cla()" -e "a.class_path=classe_path"'
else
  eval 'root  -e "std::string classe_path=\"$BASE_PATH\" " -e ".L $BASE_PATH/Class/_c/class_include.hpp" -e "#include \"test.hpp\" " -e "auto a = cla()" -e "a.class_path=classe_path"'
fi
