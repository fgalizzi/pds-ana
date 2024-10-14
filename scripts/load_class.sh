BASE_PATH="/Users/federico/pds-ana"
echo -e "\n \n Running the general class Analyzer \n \n"
eval 'root -e "#include \"const.hpp\" " -e "std::string classe_path = \"$BASE_PATH\" " -e "#include \"$BASE_PATH/Class/_c/class_include.hpp\" " -e "auto a = cla()" '
