# check for premake5
command -v premake5 >/dev/null || ./install_premake.sh

LIB_DAISY_REPO=https://github.com/electro-smith/libDaisy

if [ ! -d libDaisy ]; then
	echo "Cloning Daisy."
    git clone $LIB_DAISY_REPO
else
    echo "Already have libDaisy, update it."
    cd libDaisy && git checkout master && git fetch && git pull
    cd ..
fi

cd libDaisy
make clean && make -j7