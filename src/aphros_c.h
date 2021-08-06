// Created by Sergey Litvinov on 01.03.2021
// Copyright 2021 ETH Zurich

#ifdef __cplusplus
extern "C" {
#endif
int aphros_Main(int, const char**);
struct aphros_Parser;
struct aphros_Parser* aphros_ParserFileIni(const char*);
int aphros_ParserFin(struct aphros_Parser*);
int aphros_ParserPrintVars(struct aphros_Parser*);
int aphros_ParserStatus(struct aphros_Parser*);
int aphros_ParserGetInt(struct aphros_Parser*, const char*);
double aphros_ParserGetDouble(struct aphros_Parser*, const char*);
char* aphros_ParserGetString(struct aphros_Parser*, const char*);
double* aphros_ParserGetVect(struct aphros_Parser*, const char*, /**/ int*);
const char* aphros_GetGitRev(void);
const char* aphros_GetLogo(void);
int aphros_GetSphereOverlap(
    const double*, const double*, const double*, /**/ double*);

#ifdef __cplusplus
}
#endif
