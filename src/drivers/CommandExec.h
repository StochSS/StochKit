#ifndef COMMAND_EXEC_H
#define COMMAND_EXEC_H

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

#endif
