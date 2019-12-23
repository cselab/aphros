static int util_digits(const char* s, char* ans) {
  int i, j;

  for (i = j = 0; s[i] != '\0'; i++)
    if (s[i] == '/') j = i + 1;
  for (; s[j] != '\0'; j++)
    if (isdigit(s[j])) break;
  for (i = j; s[i] != '\0'; i++)
    if (!isdigit(s[i]))
      break;
    else
      *ans++ = s[i];
  *ans = '\0';
  return 0;
}

static int util_name(const char* p0, const char* name, char* output) {
  char dig[N], pattern[N];
  char* c;

  strncpy(pattern, p0, N - 1);
  c = strchr(pattern, '%');
  if (c == NULL) {
    fprintf(stderr, "%s: no %% in pattern '%s'\n", me, pattern);
    return 1;
  }
  *c = '\0';
  util_digits(name, dig);
  if (snprintf(output, N, "%s%s%s", pattern, dig, c + 1) < 0) {
    fprintf(stderr, "%s: snprintf failed\n", me);
    return 1;
  }
  return 0;
}

static int util_eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}
