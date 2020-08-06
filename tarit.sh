#!/usr/local/bin/zsh
# $Id: tarit.sh,v 1.7 2002/12/17 10:51:33 gsalam Exp $

release='1.0.4'
basedir=dispatch
basename=$basedir-$release



echo '---------------------------------------------'
echo 'Will create a file ../'$basename.tar.gz
echo 'Verify version number in common/run_descriptor.f90 is correct:'
grep "'Version " common/run_descriptor.f90

pushd ..

#filelist=($basedir/(*[A-Z][A-Z]|*.f|*.f90|Makefile|Makefile.???|*.fig))
#--- use the following because [A-Z] on some systems includes [a-z]...
disentdir=($basedir/disent)
disasterdir=($basedir/disaster)
dstrlibdir=($basedir/disaster/dstrlib)
commondir=($basedir/common)
exampledir=($basedir/examples)

baselist=($basedir/(CHANGELOG|README|INSTALL|MANUAL|tarit.sh))
commonlist=($commondir/(README|INSTALL|MANUAL|*.pl|*.f90|*.f))
disentlist=($disentdir/(README|INSTALL|MANUAL|*.f90|*.f|Makefile.???))
disasterlist=($disasterdir/(README|INSTALL|MANUAL|*.f90|*.f|Makefile.???))
dstrliblist=($dstrlibdir/(README|INSTALL|MANUAL|*.f|*.h|*.cc|*.f|Makefile.???))
examplelist=($exampledir/*.rdf)

tar cf - $baselist $commonlist $disentlist $disasterlist $dstrliblist $examplelist | gzip -9 > $basename.tar.gz


popd


#ls -Fd $filelist


