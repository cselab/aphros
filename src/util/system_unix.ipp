// Created by Sergey Litvinov on 08.03.2021
// Copyright 2021 ETH Zurich

#define POSIX_C_SOURCE 200809L
#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

int SystemBaseName(const char* path, char *name) {
  char *ans;
  char path0[FILENAME_MAX];
  strncpy(path0, path, sizeof path0 - 1);
  if ((ans = basename(path0)) == NULL)
    return 1;
  strcpy(name, ans);
  return 0;
}

int SystemDirName(const char *path, char *dir) {
  char path0[FILENAME_MAX];
  strncpy(path0, path, sizeof path0 - 1);
  strcpy(dir, dirname(path0));
  return 0;
}

int SystemMakeDir(const char* path, int parent) {
  char cmd[FILENAME_MAX + 200];
  int rc;

  (void)rc;
  if (parent)
    snprintf(cmd, sizeof cmd - 1, "mkdir -p -- '%s'", path);
  else
    snprintf(cmd, sizeof cmd - 1, "mkdir -- '%s'", path);
  rc = system(cmd);
  return SystemIsDir(path);
}

char* SystemRealPath(const char* path, char* resolved) {
  return realpath(path, resolved);
}

int SystemGetHostName(char* name, size_t len) {
  return gethostname(name, len);
}

int SystemHasHyperthreads(void) {
  FILE *file;
  char line[1024];
  char *search = " :\t\n";
  char *tok;
  if ((file = fopen("/proc/cpuinfo", "r")) == NULL)
    return 0;
  while (fgets(line, sizeof line - 1, file)) {
    tok = strtok(line, search);
    if (strcmp(tok, "flags") == 0) {
      while ((tok = strtok(NULL, search)) != NULL) {
        if (strcmp(tok, "ht") == 0) {
            fclose(file);
            return 1;
        }
      }
    }
  }
  fclose(file);
  return 0;
}

int SystemIsDir(const char* path) {
  struct stat info;
  return (stat(path, &info) == 0) && (info.st_mode & S_IFDIR);
}

int SystemIsFile(const char* path) {
  struct stat info;
  return (stat(path, &info) == 0) && (info.st_mode & (S_IFREG | S_IFLNK));
}

int SystemJoin(const char *a, const char *b, char *c) {
  if (a[0] == '\0')
    strcpy(c, b);
  else if (b[0] == '\0')
    strcpy(c, a);
  else if (b[0] == '/')
    strcpy(c, b);
  else if (a[strlen(a) - 1] == '/') {
    strcpy(c, a);
    strcat(c, b);
  } else {
    strcpy(c, a);
    strcat(c, "/");
    strcat(c, b);
  }
  return 0;
}

size_t SystemGetMem(void) {
  FILE *file;
  char line[1024];
  char *search = " :\t\n";
  char *tok;
  char *end;
  size_t ans;
  if ((file = fopen("/proc/self/status", "r")) == NULL)
    return 0;
  while (fgets(line, sizeof line - 1, file)) {
    tok = strtok(line, search);
    if (strcmp(tok, "VmRSS") == 0) {
      if ((tok = strtok(NULL, search)) == NULL) {
        fclose(file);
        return 0;
      }
      ans = strtol(tok, &end, 10);
      fclose(file);
      return ans << 10;
    }
  }
  fclose(file);
  return 0;
}

int SystemSplitExt(const char *a, char *base, char *ext) {
  int i;
  int j;
  int k;
  for (i = 0, j = -1; a[i] != '\0'; i++)
    if (a[i] == '/') j = i;
  for (i = j + 1, k = 0; a[i] != '\0'; i++)
    if (a[i] == '.') k = i;
  if (k == 0 || (j != -1 && k == j + 1)) {
    strcpy(base, a);
    ext[0] = '\0';
  } else {
    strncpy(base, a, k);
    base[k] = '\0';
    strcpy(ext, &a[k]);
  }
  return 0;
}
