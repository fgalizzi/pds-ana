C++ class to analyze DUNE PDS data

Call it with "c" from everywhere on your terminal

If you want to modify it, copy and paste it and the folder "_ c" wherever
you want; then, call it with "cx"



.zshrc aliases

alias c="source ~/.load_class.sh"
  eval 'root -e "#include \"/Users/federico/PhD/Class/classe.hpp\" " -e "auto a = cla()" '
alias cx="source ~/.load_class_here.sh"
  eval 'root -e "#include \"classe.hpp\" " -e "auto a = cla()" '
