function [lf] = gradllfCovarpara_nogradient(parametervectorPre,ntau, nsample, nmixtures, nslice, yt,X)

[lf,g] = gradllfCovarpara(parametervectorPre,ntau, nsample, nmixtures, nslice, yt,X)
