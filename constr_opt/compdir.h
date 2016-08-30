
#ifndef _COMPDIR_H_
#define _COMPDIR_H_


struct CompDirRes
{
	EVector SD;
	SearchDir dirType;
};

CompDirRes CompDir(const EMatrix* Z, const EMatrix* H, const EVector* gf, int numVars, const EVector* f, float eps);

#endif
