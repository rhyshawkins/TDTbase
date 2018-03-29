#ifndef slog_h
#define slog_h

#include <stdio.h>
#include <stdarg.h>

enum {
  SLOG_FLAGS_NONE = 0,
  SLOG_FLAGS_CLEAR = 1
};

enum {
  SLOG_ERROR,
  SLOG_WARNING,
  SLOG_INFO,
  SLOG_DEBUG
};

/*
 * By default, logs go to stderr, this redirects to a file which is optionally overwritten.
 */
int slog_set_output_file(const char *filename,
			 int flags);

#define ERROR(fmt, ...) slog(SLOG_ERROR, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define WARNING(fmt, ...) slog(SLOG_WARNING, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define INFO(fmt, ...) slog(SLOG_INFO, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define DEBUG(fmt, ...) slog(SLOG_WARNING, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

int slog(int level,
	 const char *source_file,
	 const char *function,
	 int lineno,
	 const char *msg,
	 ...);

int vslog(int level,
	  const char *source_file,
	  const char *function,
	  int lineno,
	  const char *msg,
	  va_list ap);

#define INFO_LARGE_START(fmt, ...) slog_large_message_start(SLOG_INFO, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define INFO_LARGE_WRITE(fmt, ...) slog_large_message_write(fmt, ##__VA_ARGS__)
#define INFO_LARGE_NEWLINE() slog_large_message_newline()
#define INFO_LARGE_END() slog_large_message_end()

int slog_large_message_start(int level,
			     const char *source_file,
			     const char *function,
			     int lineno,
			     const char *msg,
			     ...);

int slog_large_message_write(const char *msg, ...);

int slog_large_message_newline(void);

int slog_large_message_end(void);

#define SLOGSTART() slogstart(__FILE__, __FUNCTION__, __LINE__)
#define SLOGEND() slogend()

FILE *slogstart(const char *source_file,
		 const char *function,
		 int lineno);
void slogend();


#endif /* log_h */
