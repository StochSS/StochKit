#include "CommandPass.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fcntl.h>
#include <sys/wait.h>

namespace STOCHKIT
{

static inline void checkError(const long int retCode, const char *const msg)
{
   if(retCode != 0)
   {
      std::cerr << "error: " << msg << std::endl;
      perror("");
      exit(1);
   }
}

static inline void checkError2(const long int retCode, const char *const msg)
{
   if(retCode == -1)
   {
      std::cerr << "error: " << msg << std::endl;
      perror("");
      exit(1);
   }
}

static inline void checkError3(const long int retCode, const char *const msg)
{
   if(retCode < 0)
   {
      std::cerr << "error: " << msg << std::endl;
      exit(1);
   }
}

CommandPass::CommandPass()
{
   processActive = false;
   start();
   //std::cout << "Proc init" << std::endl;
}

void CommandPass::start(void)
{
   int cmd_fd[2], stat_fd[2];
   checkError(pipe(cmd_fd), "pipe error");
   checkError(pipe(stat_fd), "pipe error");
   proc = fork();
   checkError2(proc, "fork error");
   if(proc == 0)
   {
      checkError(close(cmd_fd[1]), "close error");
      checkError(close(stat_fd[0]), "close error");
      char buf1[100], buf2[100];
      char bufPath[10020];
      const char *path_ptr;
      path_ptr = getenv("STOCHKIT_HOME");
      if(path_ptr == NULL)
      {
         std::cerr << "STOCHKIT_HOME does not exist" << std::endl;
         exit(1);
      }
      memset(bufPath, 0, 10020);
      strncpy(bufPath, path_ptr, 10000);
      strncat(bufPath, "/bin/CommandExec", 17);
      checkError3(snprintf(buf1, 100, "%d", cmd_fd[0]), "snprintf error");
      checkError3(snprintf(buf2, 100, "%d", stat_fd[1]), "snprintf error");
      //std::cout << bufPath << ";" << buf1 << ";" << buf2 << std::endl;
      checkError2(execl(bufPath, "CommandExec", buf1, buf2, NULL), "exec error");
   }
   else
   {
      processActive = true;
      pipeFD = cmd_fd[1];
      statusFD = stat_fd[0];
      checkError2(fcntl(cmd_fd[1], F_SETFD, FD_CLOEXEC), "fcntl error");
      checkError2(fcntl(stat_fd[0], F_SETFD, FD_CLOEXEC), "fcntl error");
      checkError(close(cmd_fd[0]), "close error");
      checkError(close(stat_fd[1]), "close error");
   }
}

static inline void readFully(const int fildes, void *const buf, const size_t nbyte)
{
   size_t total = 0;
   ssize_t ret;
   while(total < nbyte)
   {
      ret = read(fildes, reinterpret_cast<char *>(buf) + total, nbyte - total);
      checkError2(ret, "read error");
      total += static_cast<size_t>(ret);
   }
}

int CommandPass::execute(const char *const cmd)
{
   size_t counter;
   int retVal;
   counter = strlen(cmd) + 1;
   checkError2(write(pipeFD, &counter, sizeof(counter) ), "write error");
   checkError2(write(pipeFD, cmd, counter), "write error");
   readFully(statusFD, &retVal, sizeof(retVal) );
   return retVal;
}

void CommandPass::stop(void)
{
   if(processActive == true)
   {
      int childStatus;
      size_t counter;
      counter = 2000000000;
      checkError2(write(pipeFD, &counter, sizeof(counter) ), "write error");
      checkError2(waitpid(proc, &childStatus, 0), "waitpid error");
      if(!WIFEXITED(childStatus) || (WIFEXITED(childStatus) && WEXITSTATUS(childStatus) != 0) )
      {
         std::cerr << "Child process error" << std::endl;
         exit(1);
      }
      checkError(close(pipeFD), "close error");
      checkError(close(statusFD), "close error");
      processActive = false;
   }
}

CommandPass::~CommandPass()
{
   stop();
}

}
