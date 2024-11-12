BASE_PATH="YOUR_PATH_TO/pds-ana"
echo -e "\n \n Running the general class Analyzer \n \n"
eval 'root  -e "#include \"const.hpp\" " -e "std::string classe_path=\"$BASE_PATH\" " -e "#include \"$BASE_PATH/Class/_c/class_include.hpp\" " -e "auto a = cla()" -e "a.class_path=classe_path"'
