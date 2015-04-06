#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string>

#include "CommandExec.h"

int main(const int argc, const char *const argv[])
{
   if(argc < 3)
   {
      std::cerr << "Need more arguments" << std::endl;
      exit(1);
   }
   int cmd_fd, stat_fd;
   int retVal;
   size_t counter;
   std::string cmd;
   
   cmd_fd = atoi(argv[1]);
   stat_fd = atoi(argv[2]);
   checkError2(fcntl(cmd_fd, F_SETFD, FD_CLOEXEC), "fcntl error");
   checkError2(fcntl(stat_fd, F_SETFD, FD_CLOEXEC), "fcntl error");

   for(;;)
   {
      cmd.clear();
      readFully(cmd_fd, &counter, sizeof(counter) );
      if(counter >= 2000000000)
      {
         break;
      }
      cmd.resize(counter);
      readFully(cmd_fd, &cmd.at(0), counter);
      retVal = system(cmd.c_str() );
      checkError2(write(stat_fd, &retVal, sizeof(retVal) ), "write error");
   }
   checkError(close(cmd_fd), "close error");
   checkError(close(stat_fd), "close error");
   return 0;
}
