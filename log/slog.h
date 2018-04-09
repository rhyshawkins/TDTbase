//
//    Simply Logging Library  
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
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
