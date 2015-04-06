#ifndef COMMAND_PASS_H
#define COMMAND_PASS_H

#include <unistd.h>

namespace STOCHKIT
{

class CommandPass
{
private:
   bool processActive;
   pid_t proc;
   int pipeFD, statusFD;

public:
   
   CommandPass();
   void start(void);
   int execute(const char *const cmd);
   void stop(void);
   ~CommandPass();
};

}

#endif
