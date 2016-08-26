
#ifndef _COMPDIR_H_
#define _COMPDIR_H_

enum SearchDir
{
	kNewton,
	kRandom,
	kEigenvec,
	kStpDesc,
};

struct CompDirRes
{
	EVector SD;
	int dirType;
};

CompDirRes CompDir(const EMatrix* Z, const EMatrix* H, const EVector* gf, int numVars, const EVector* f, float eps);

#endif
