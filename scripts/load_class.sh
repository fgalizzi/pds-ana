echo -e "\n \n Running the general class Analyzer \n \n"
eval 'root -e "#include \"const.hpp\" " -e "std::string classe_path = \"/Users/federico/pds-ana/Class\"" -e "#include \"/Users/federico/pds-ana/Class/classe.hpp\" " -e "auto a = cla()" '
