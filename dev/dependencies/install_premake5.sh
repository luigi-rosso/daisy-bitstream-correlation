#!/bin/sh
echo Installing Premake5
wget https://github.com/premake/premake-core/releases/download/v5.0.0-alpha15/premake-5.0.0-alpha15-linux.tar.gz
tar -xvf premake-5.0.0-alpha15-linux.tar.gz
mkdir bin
cp premake5 bin/premake5
mv premake5 /usr/local/bin