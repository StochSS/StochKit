/*
*  shell_tau_leaping.cpp
*  call .\bin\tau_leaping
*/

#include "boost_headers.h"

int main(int ac, char* av[])
{
	boost::filesystem::path currentPath=boost::filesystem::system_complete(av[0]).parent_path(); 
	std::string command=currentPath.string();
	command+="\\bin\\tau_leaping_hub";

	//add quote for the executable
	command="\""+command+"\"";

	for(int i=1; i<ac; i++)
	{
		if((std::string)av[i]=="-m" || (std::string)av[i]=="--out-dir" || (std::string)av[i]=="--model")
		{
			command+=" "+std::string(av[i]);
			command+=" \""+std::string(av[i+1])+"\"";
			i++;
		}
		else
			command+=" "+std::string(av[i]);
	}  
	//add quote for the command
	command="\""+command+"\""; 
	system(command.c_str());

	return 0;
}
