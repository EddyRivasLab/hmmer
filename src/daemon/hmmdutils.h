extern void p7_openlog(const char *ident, int option, int facility);
extern void p7_syslog(int priority, const char *format, ...);
extern void p7_closelog(void);
